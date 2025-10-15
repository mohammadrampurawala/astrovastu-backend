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
from dateutil import parser as dateutil_parser
import pytz
from timezonefinder import TimezoneFinder
from geopy.geocoders import Nominatim
from fastapi.responses import JSONResponse
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

# -----------------------
# Helper & Request Model
# -----------------------
class BirthData(BaseModel):
    date: str
    time: str
    place: str
    lat: float | None = None
    lon: float | None = None
    timezone: str | None = None
    sidereal: bool | None = False

def normalize_angle(deg):
    """Normalize angle to 0..360 float."""
    try:
        d = float(deg) % 360.0
        if d < 0:
            d += 360.0
        return d
    except Exception:
        return None

def _extract_from_res(res):
    """
    Accepts res which might be:
     - a float/number
     - a tuple/list with elements (longitude, latitude, distance, speed_long, ...)
    Returns dict with keys: longitude, latitude, speed_long (values may be None).
    """
    out = {"longitude": None, "latitude": None, "speed_long": None}
    try:
        # sequence-like
        if isinstance(res, (list, tuple)):
            if len(res) >= 1:
                out["longitude"] = float(res[0])
            if len(res) >= 2:
                out["latitude"] = float(res[1])
            if len(res) >= 4:
                out["speed_long"] = float(res[3])
        else:
            # scalar case (just a number)
            out["longitude"] = float(res)
    except Exception as e:
        # log for debugging; Render will show prints
        print(f"Warning: unable to parse ephemeris result {res!r} -> {e}")
    return out

def _parse_datetime_to_utc_jd(date_str: str, time_str: str, tz_name: str | None):
    """
    Parse local date/time to UTC datetime and return Swiss Ephemeris Julian Day (UT).
    """
    # Combine date/time and parse
    dt_local = dateutil_parser.parse(f"{date_str} {time_str}")
    # If user supplied timezone, use it; else assume naive -> treat as local system tz?
    # Prefer explicit timezone input; if missing, assume UTC for safety (or you can default to Asia/Kolkata)
    if tz_name:
        try:
            tz = pytz.timezone(tz_name)
            dt_local = tz.localize(dt_local) if dt_local.tzinfo is None else dt_local.astimezone(tz)
        except Exception as e:
            # fallback: assume naive datetimes are local system tz
            print(f"Warning: timezone parse failed ({tz_name}) -> {e}")
    # Ensure dt is timezone-aware and convert to UTC
    if dt_local.tzinfo is None:
        # treat naive as UTC to avoid ambiguity; if you prefer treat as local, change here
        dt_utc = dt_local.replace(tzinfo=pytz.UTC)
    else:
        dt_utc = dt_local.astimezone(pytz.UTC)
    # compute decimal hour
    h = dt_utc.hour + dt_utc.minute / 60.0 + dt_utc.second / 3600.0 + dt_utc.microsecond / 3_600_000_000.0
    # use swe.julday
    jd_ut = swe.julday(dt_utc.year, dt_utc.month, dt_utc.day, h)
    return dt_utc, jd_ut

def _parse_place_to_latlon(place_str: str):
    """
    Accept a place string in the form "lat,lon" and return (lat, lon).
    If not parseable, return (None, None) — caller may try geocoding if desired.
    """
    try:
        parts = [p.strip() for p in place_str.split(",")]
        if len(parts) >= 2:
            lat = float(parts[0])
            lon = float(parts[1])
            # sanity check ranges
            if -90 <= lat <= 90 and -180 <= lon <= 180:
                return lat, lon
    except Exception:
        pass
    return None, None

# -----------------------
# compute_chart endpoint
# -----------------------
@app.post("/compute_chart")
def compute_chart(payload: BirthData):
    """
    Compute planetary longitudes, house cusps and ascendant.
    Returns an object containing 'planets', 'houses', 'ascendant', and input/utc_birth.
    This implementation is defensive and will return helpful error JSON if required fields cannot be computed.
    """
    data = payload.dict()
    # 1) parse lat/lon from explicit fields or place string "lat,lon"
    lat = payload.lat
    lon = payload.lon
    if lat is None or lon is None:
        latlon = _parse_place_to_latlon(payload.place)
        if latlon != (None, None):
            lat, lon = latlon

    if lat is None or lon is None:
        # We deliberately require lat/lon or a lat,lon place string. If you want geocoding, add geopy.
        return JSONResponse(status_code=422, content={
            "error": "missing_latlon",
            "message": "Please provide latitude and longitude either via 'lat' and 'lon' fields or as 'place' in the form 'lat,lon'."
        })

    # 2) parse date/time -> UTC and JD
    try:
        dt_utc, jd_ut = _parse_datetime_to_utc_jd(payload.date, payload.time, payload.timezone)
    except Exception as e:
        return JSONResponse(status_code=422, content={
            "error": "invalid_datetime",
            "message": str(e)
        })

    # 3) Build planets
    # Use common SWEP constants
    planet_ids = {
        "Sun": swe.SUN,
        "Moon": swe.MOON,
        "Mercury": swe.MERCURY,
        "Venus": swe.VENUS,
        "Mars": swe.MARS,
        "Jupiter": swe.JUPITER,
        "Saturn": swe.SATURN,
        # use mean lunar node for Rahu (name as 'Rahu') and Ketu as opposite point
        "Rahu": swe.MEAN_NODE,
        "Ketu": swe.TRUE_NODE  # if you prefer Ketu = Rahu+180 adjust later
    }

    planets_out = {}
    for name, pid in planet_ids.items():
        try:
            # call Swiss Ephemeris; calc_ut may return tuple/list or scalar depending on flags/version
            res = swe.calc_ut(jd_ut, pid)
        except Exception as e:
            print(f"Ephemeris call failed for {name}: {e}")
            res = None

        vals = _extract_from_res(res)
        if vals["longitude"] is None:
            # log and continue (do not crash)
            print(f"compute_chart: missing longitude for {name}; raw res={res!r}")
            planets_out[name] = {
                "longitude": None,
                "latitude": vals.get("latitude"),
                "speed_long": vals.get("speed_long")
            }
        else:
            lon_deg = normalize_angle(vals["longitude"])
            planets_out[name] = {
                "longitude": lon_deg,
                "latitude": vals.get("latitude"),
                "speed_long": vals.get("speed_long")
            }

    # 4) Houses & Ascendant
    # Swiss Ephemeris: swe.houses(jd_ut, lat, lon) -> (cusps, ascmc)
    try:
        cusps, ascmc = swe.houses(jd_ut, lat, lon)
        # cusps is an array-like with 13 entries (1..12)
        houses = {}
        # ensure numeric keys "1".."12"
        for i in range(1, 13):
            try:
                houses[str(i)] = normalize_angle(cusps[i])
            except Exception:
                houses[str(i)] = None
        # ascendant value usually ascmc[0]
        asc_value = None
        try:
            asc_value = normalize_angle(ascmc[0])
        except Exception:
            asc_value = None
    except Exception as e:
        print(f"House calculation failed: {e}")
        houses = {str(i): None for i in range(1, 13)}
        asc_value = None

    # 5) Validate minimal outputs (planets dict must have Jupiter, houses must include 5th)
    missing = []
    if "Jupiter" not in planets_out or planets_out["Jupiter"]["longitude"] is None:
        missing.append("planet:Jupiter")
    if houses.get("5") is None:
        missing.append("house:5")
    if asc_value is None:
        missing.append("ascendant")

    if missing:
        # return helpful diagnostic JSON for GPT and debugging
        return JSONResponse(status_code=500, content={
            "error": "incomplete_chart",
            "missing_fields": missing,
            "input": data,
            "planets_sample": {k: v for k, v in list(planets_out.items())[:5]},
            "houses_sample": {k: houses[k] for k in list(houses.keys())[:6]}
        })

    # 6) Build and return response
    response = {
        "input": data,
        "utc_birth": dt_utc.isoformat(),
        "planets": planets_out,
        "ascendant": asc_value,
        "houses": houses
    }
    return response


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




