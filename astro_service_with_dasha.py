# astro_service_with_dasha.py
"""
FastAPI service that computes natal chart data and Vimshottari mahadasha + antardasha timelines
using Swiss Ephemeris (pyswisseph).

Endpoints:
 - POST /compute_chart  -> returns planets, ascendant, houses, utc_birth, jd_ut
 - POST /compute_dasha  -> returns moon longitude, nakshatra, mahadasha sequence and nested antardashas

Install dependencies:
 pip install pyswisseph fastapi uvicorn geopy timezonefinder python-dateutil

Run:
 uvicorn astro_service_with_dasha:app --port 8000 --reload
"""
import swisseph as swe
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from datetime import datetime, timezone, timedelta
from dateutil import tz
from timezonefinder import TimezoneFinder
from geopy.geocoders import Nominatim
import math

app = FastAPI(title="Astro Service with Dasha")

# Set ephemeris path if you have local ephe files, otherwise default built-ins
swe.set_ephe_path('.')

# Planet constants map
PLANETS = {
    'Sun': swe.SUN,
    'Moon': swe.MOON,
    'Mercury': swe.MERCURY,
    'Venus': swe.VENUS,
    'Mars': swe.MARS,
    'Jupiter': swe.JUPITER,
    'Saturn': swe.SATURN,
    'Rahu': swe.MEAN_NODE,  # use mean node, treat Ketu separately
}

# Vimshottari standard order and durations (years)
VIMSHOTTARI_ORDER = [
    "Ketu", "Venus", "Sun", "Moon", "Mars",
    "Rahu", "Jupiter", "Saturn", "Mercury"
]
VIMSHOTTARI_YEARS = {
    "Ketu": 7.0, "Venus": 20.0, "Sun": 6.0, "Moon": 10.0,
    "Mars": 7.0, "Rahu": 18.0, "Jupiter": 16.0, "Saturn": 19.0, "Mercury": 17.0
}
TOTAL_CYCLE_YEARS = sum(VIMSHOTTARI_YEARS.values())  # should be 120

NAK_LEN = 360.0 / 27.0  # 13.333333333333334

# Mapping nakshatra index (0..26) to starting planet index in VIMSHOTTARI_ORDER
NAK_TO_START_PLANET_INDEX = [
    0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8
]

geolocator = Nominatim(user_agent="astrovastu_pro")
tzfinder = TimezoneFinder()

class BirthData(BaseModel):
    date: str  # YYYY-MM-DD
    time: str  # HH:MM (24h)
    place: str  # "City, Country" or "lat,lon"
    lat: float = None
    lon: float = None
    timezone: str = None  # optional IANA
    sidereal: bool = False  # default tropical

def normalize_angle(angle):
    """
    Normalize an angle to 0..360 degrees.
    Accepts a float/int OR a sequence (tuple/list/np.array) where the first element is the angle.
    """
    # If angle is a sequence (tuple/list/ndarray...), try to extract the first element
    try:
        # Protect against strings (they are sequences but not what we want)
        if not isinstance(angle, (str, bytes)) and hasattr(angle, "__len__"):
            # get first element (works for tuple/list/numpy array)
            angle_val = angle[0]
        else:
            angle_val = angle
        a = float(angle_val)
    except Exception:
        # Re-raise with helpful message for debugging
        raise TypeError(f"normalize_angle: could not convert angle {angle!r} to float")

    a = a % 360.0
    if a < 0:
        a += 360.0
    return a

def parse_place(place, lat, lon):
    if lat is not None and lon is not None:
        return lat, lon
    try:
        if ',' in place:
            parts = place.split(',')
            latf = float(parts[0].strip()); lonf = float(parts[1].strip())
            return latf, lonf
    except:
        pass
    loc = geolocator.geocode(place, timeout=10)
    if not loc:
        raise ValueError("Could not geocode place. Provide lat,lon or clearer place.")
    return loc.latitude, loc.longitude

def local_to_utc_datetime(date_str, time_str, tz_name, lat=None, lon=None):
    """Return a timezone-aware UTC datetime for the provided local date/time & tz name.
       tz_name if None will be inferred from lat/lon."""
    naive = datetime.strptime(f"{date_str} {time_str}", "%Y-%m-%d %H:%M")
    if not tz_name:
        if lat is None or lon is None:
            raise ValueError("Timezone unknown and lat/lon not provided.")
        tz_name = tzfinder.timezone_at(lat=lat, lng=lon)
        if not tz_name:
            raise ValueError("Could not infer timezone.")
    local_tz = tz.gettz(tz_name)
    if local_tz is None:
        raise ValueError("Invalid timezone provided.")
    local_dt = naive.replace(tzinfo=local_tz)
    utc_dt = local_dt.astimezone(timezone.utc)
    return utc_dt

def jd_from_utc_dt(utc_dt):
    """Compute Julian day (UT) from a timezone-aware UTC datetime."""
    year = utc_dt.year; month = utc_dt.month; day = utc_dt.day
    hour = utc_dt.hour + utc_dt.minute/60.0 + utc_dt.second/3600.0
    return swe.julday(year, month, day, hour)

@app.post("/compute_chart")
def compute_chart(data: BirthData):
    lat, lon = parse_place(data.place, data.lat, data.lon)
    tz_name = data.timezone or tzfinder.timezone_at(lat=lat, lng=lon)
    try:
        utc_dt = local_to_utc_datetime(data.date, data.time, tz_name, lat, lon)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

    jd_ut = jd_from_utc_dt(utc_dt)

    # set sidereal mode if requested (Lahiri)
    if data.sidereal:
        swe.set_sid_mode(swe.SIDM_LAHIRI)
    else:
        swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY)

    # compute planets
    planets_out = {}
    for name, pid in PLANETS.items():
        res = swe.calc_ut(jd_ut, pid)
        lon_deg = normalize_angle(res[0])
        planets_out[name] = {'longitude': lon_deg, 'latitude': res[1], 'speed_long': res[3]}

    # compute Ketu as Rahu + 180
    rahu_lon = planets_out['Rahu']['longitude']
    ketu_lon = normalize_angle(rahu_lon + 180.0)
    planets_out['Ketu'] = {'longitude': ketu_lon, 'latitude': None, 'speed_long': None}

    # houses & ascendant
    hsys = 'P'
    cusps, ascmc = swe.houses(jd_ut, lat, lon, hsys)
    asc = normalize_angle(ascmc[0])
    mc = normalize_angle(ascmc[1])
    houses = {f"house_{i+1}": cusps[i] for i in range(12)}

    return {
        'input': {'date': data.date, 'time': data.time, 'place': data.place, 'lat': lat, 'lon': lon, 'timezone': tz_name, 'sidereal': data.sidereal},
        'utc_birth': utc_dt.isoformat(),
        'jd_ut': jd_ut,
        'planets': planets_out,
        'ascendant': asc,
        'mc': mc,
        'houses': houses
    }

# ---- Dasha / Nakshatra utilities ----
def moon_to_nakshatra_index(moon_longitude_deg):
    moon = normalize_angle(moon_longitude_deg)
    nak_index_float = moon / NAK_LEN
    nak_index = int(math.floor(nak_index_float))  # 0..26
    frac_in_nak = nak_index_float - nak_index
    return nak_index, frac_in_nak

def build_mahadasha_sequence(moon_lon_deg, birth_utc_dt):
    nak_index, frac = moon_to_nakshatra_index(moon_lon_deg)
    start_idx = NAK_TO_START_PLANET_INDEX[nak_index]
    ordered_planets = []
    for i in range(len(VIMSHOTTARI_ORDER)):
        ordered_planets.append(VIMSHOTTARI_ORDER[(start_idx + i) % len(VIMSHOTTARI_ORDER)])

    seq = []
    running_start = birth_utc_dt
    for idx, planet in enumerate(ordered_planets):
        full_years = VIMSHOTTARI_YEARS[planet]
        if idx == 0:
            years = full_years * (1.0 - frac)
        else:
            years = full_years
        days = years * 365.2425
        running_end = running_start + timedelta(days=days)
        seq.append({
            'planet': planet,
            'start_utc': running_start.isoformat(),
            'end_utc': running_end.isoformat(),
            'duration_years': years
        })
        running_start = running_end
    return {
        'nakshatra_index': nak_index,
        'nakshatra_fraction': frac,
        'mahadasha_sequence': seq
    }

def build_antardashas_for_mahadasha(maha_planet, maha_start_dt, maha_years):
    """Return a list of antardashas (subperiods) within the mahadasha period.
       Order of antardashas starts with maha_planet and follows VIMSHOTTARI_ORDER.
    """
    # find start index in global order
    try:
        start_index = VIMSHOTTARI_ORDER.index(maha_planet)
    except ValueError:
        start_index = 0
    antardashas = []
    maha_days = maha_years * 365.2425
    running_start = maha_start_dt
    for i in range(len(VIMSHOTTARI_ORDER)):
        sub_index = (start_index + i) % len(VIMSHOTTARI_ORDER)
        subplanet = VIMSHOTTARI_ORDER[sub_index]
        # subperiod length = maha_years * (subplanet_years / TOTAL_CYCLE_YEARS)
        sub_years = maha_years * (VIMSHOTTARI_YEARS[subplanet] / TOTAL_CYCLE_YEARS)
        sub_days = sub_years * 365.2425
        running_end = running_start + timedelta(days=sub_days)
        antardashas.append({
            'planet': subplanet,
            'start_utc': running_start.isoformat(),
            'end_utc': running_end.isoformat(),
            'duration_years': sub_years
        })
        running_start = running_end
    return antardashas

@app.post("/compute_dasha")
def compute_dasha(data: BirthData):
    # reuse compute_chart minimal steps to get moon lon and birth utc dt
    lat, lon = parse_place(data.place, data.lat, data.lon)
    tz_name = data.timezone or tzfinder.timezone_at(lat=lat, lng=lon)
    try:
        utc_dt = local_to_utc_datetime(data.date, data.time, tz_name, lat, lon)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    jd_ut = jd_from_utc_dt(utc_dt)

    # sidereal?
    if data.sidereal:
        swe.set_sid_mode(swe.SIDM_LAHIRI)
    else:
        swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY)

    # get Moon longitude (UT)
    moon_res = swe.calc_ut(jd_ut, swe.MOON)
    moon_lon = normalize_angle(moon_res[0])

    maha_info = build_mahadasha_sequence(moon_lon, utc_dt)

    # Now compute antardashas for each mahadasha entry
    # We'll also include a flattened antardasha list keyed by mahadasha
    for i, maha in enumerate(maha_info['mahadasha_sequence']):
        # parse start as datetime object
        start_dt = datetime.fromisoformat(maha['start_utc'])
        duration_years = maha['duration_years']
        maha['antardashas'] = build_antardashas_for_mahadasha(maha['planet'], start_dt, duration_years)

    result = {
        'input': {'date': data.date, 'time': data.time, 'place': data.place, 'lat': lat, 'lon': lon, 'timezone': tz_name, 'sidereal': data.sidereal},
        'utc_birth': utc_dt.isoformat(),
        'moon_longitude': moon_lon,
        'nakshatra_index': maha_info['nakshatra_index'],
        'nakshatra_fraction': maha_info['nakshatra_fraction'],
        'mahadasha_sequence': maha_info['mahadasha_sequence']
    }
    return result

# To run:
# uvicorn astro_service_with_dasha:app --port 8000 --reload
