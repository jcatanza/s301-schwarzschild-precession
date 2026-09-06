"""
Solar-conjunction geometry of Sgr A* around S301's first campaign
periapsis (k.NEXT_PERIAPSIS_YR, late October 2031).

campaign.py encodes Paranal's Sgr A* observing season as a fixed
fractional-year window (campaign._in_visibility_season) and reports that
periapsis 1 falls just outside it, so the dense window for that passage
is anchored to the last visible date instead. This script puts real
numbers behind that statement:

  - the actual solar elongation of Sgr A* (RA 17h45m40.04s,
    Dec -29d00m28.1s), daily from 2031-08-01 to 2032-04-30, via
    astropy.coordinates.get_sun -- the dates the elongation drops below
    85 deg (roughly where a ~2-3 hour VLTI observing window is lost) and
    below 45 deg (roughly unobservable), its minimum, and the elongation
    at each of the two campaign periapsis passages;
  - the campaign's own season filter, scanned daily: the last visible
    date before periapsis 1, the first visible date after, and the gap
    in days that the filter actually imposes;
  - the calendar drift, P mod 1 yr, that moves the second passage
    (~2040.49, late June) safely into the season.

ERFA emits "dubious year" warnings for UTC dates beyond its leap-second
table; they are harmless for a degree-level elongation and suppressed.

Results go to results/solar_conjunction.json (dates as 'YYYY-MM-DD').
"""

import warnings

import numpy as np
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time, TimeDelta

import campaign
import constants as k
import results_io

warnings.filterwarnings("ignore", message=".*dubious year.*")

SGRA = SkyCoord("17h45m40.04s", "-29d00m28.1s", frame="icrs")
SCAN_START = "2031-08-01"
SCAN_END = "2032-04-30"
SEASON_SCAN_HALF_DAYS = 200


def daily_grid(start, end):
    """Daily astropy Time grid (00:00 UTC) from `start` to `end` inclusive."""
    t_start, t_end = Time(start), Time(end)
    n_days = int(round((t_end - t_start).jd)) + 1
    return t_start + TimeDelta(np.arange(n_days), format="jd")


def elongation_deg(times):
    """Sun -> Sgr A* angular separation (deg) as seen from Earth at
    `times`. Sgr A* is transformed into the Sun's geocentric frame (a
    source at infinity: aberration only), so the mismatch of frame
    origins is intentional and its warning suppressed."""
    sep = get_sun(times).separation(SGRA, origin_mismatch="ignore")
    return np.atleast_1d(sep.deg)


def date_str(time):
    """'YYYY-MM-DD' for a scalar Time."""
    return time.iso[:10]


def threshold_window(times, elong, threshold_deg):
    """First and last scan dates on which the elongation is below
    threshold_deg."""
    idx = np.flatnonzero(elong < threshold_deg)
    return date_str(times[idx[0]]), date_str(times[idx[-1]])


def season_gap(peri_time):
    """Scan the campaign's own visibility filter day by day around a
    periapsis time: (last visible date before, first visible date after)."""
    days = np.arange(-SEASON_SCAN_HALF_DAYS, SEASON_SCAN_HALF_DAYS + 1)
    times = Time(date_str(peri_time)) + TimeDelta(days, format="jd")
    # pylint: disable=protected-access
    # The season filter is deliberately private to campaign.py (the grid
    # builder is its only intended caller); this script exists to audit
    # exactly that filter, so it must call the same function.
    visible = campaign._in_visibility_season(times.decimalyear)
    # pylint: enable=protected-access
    before = times[visible & (times < peri_time)]
    after = times[visible & (times > peri_time)]
    return before[-1], after[0]


def main():
    """Daily elongation scan, season-filter gap, and results."""
    times = daily_grid(SCAN_START, SCAN_END)
    elong = elongation_deg(times)
    j_min = int(np.argmin(elong))
    peri_one = Time(k.NEXT_PERIAPSIS_YR, format="decimalyear")
    peri_two = Time(k.NEXT_PERIAPSIS_YR + k.TRUTH["P_yr"], format="decimalyear")
    peri_one_elong = float(elongation_deg(peri_one)[0])
    peri_two_elong = float(elongation_deg(peri_two)[0])

    results = {
        "elong_min_deg": (float(elong[j_min]), ".1f"),
        "elong_min_date": date_str(times[j_min]),
        "peri_one_date": date_str(peri_one),
        "peri_one_elong_deg": (peri_one_elong, ".0f"),
        "peri_two_elong_deg": (peri_two_elong, ".0f"),
        "calendar_drift_yr": (k.TRUTH["P_yr"] % 1, ".2f"),
    }
    print(f"Sgr A* solar elongation, daily {SCAN_START} to {SCAN_END}: minimum {elong[j_min]:.1f} deg "
          f"on {results['elong_min_date']}")
    for threshold in (45.0, 85.0):
        start, end = threshold_window(times, elong, threshold)
        results[f"elong_below{threshold:.0f}_start"] = start
        results[f"elong_below{threshold:.0f}_end"] = end
        print(f"  elongation < {threshold:.0f} deg: {start} to {end}")
    print(f"Periapsis 1 at {peri_one.iso} (decimal year {k.NEXT_PERIAPSIS_YR:.3f}): "
          f"elongation {peri_one_elong:.0f} deg")
    print(f"Periapsis 2 at {peri_two.iso}: elongation {peri_two_elong:.0f} deg "
          f"(calendar drift P mod 1 yr = {k.TRUTH['P_yr'] % 1:.2f} yr)")

    last_before, first_after = season_gap(peri_one)
    gap_before = float((peri_one - last_before).jd)
    gap_after = float((first_after - peri_one).jd)
    gap_total = float((first_after - last_before).jd)
    print(f"campaign._in_visibility_season: last visible date before periapsis 1 "
          f"{date_str(last_before)} ({gap_before:.0f} d before), first visible date after "
          f"{date_str(first_after)} ({gap_after:.0f} d after); total gap {gap_total:.0f} d")
    results.update({
        "gap_days_before_peri": (gap_before, ".0f"),
        "gap_days_after_peri": (gap_after, ".0f"),
        "gap_total_days": (gap_total, ".0f"),
    })
    results_io.write_results("solar_conjunction", results)


if __name__ == "__main__":
    main()
