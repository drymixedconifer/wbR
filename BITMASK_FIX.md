# v6 monthly bitmask correction

* v6 has 14 monthly variables, including runoff at bit 8192, so accepted masks are 0–16383 (previous validator erroneously stopped at 8191).
* Runoff column cursor advances consistently.
* Zonal-only simulations request only needed monthly fields rather than all 14 unconditionally, reducing memory demand.
* No water-balance equations or PRISM input processing changed.
* Source inspected; R executable unavailable in build environment.
