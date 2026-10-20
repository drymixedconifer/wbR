# Seasonal snow initialization correction

The first baseline year previously created an all-NA `snow`, `pack`, or `snow_peak` band because the spin-up checkpoint held December soil and SWE but did not hold the preceding September--December monthly SWE history.

The first annual sequence now uses one additional PRISM-normal simulation cycle if the preceding autumn is absent. It carries that cycle's December soil/snow state forward along with the September--December SWE sum and peak. This preserves the nine-month September--May definition for the first output year, and also supports analysis initialized from equilibrium rather than baseline end. Later years continue to derive their preceding autumn from the actual previous simulated year.

No snow equations were changed. One additional normal year may slightly adjust the initial soil/snow state; compare first-year outputs when evaluating backward compatibility. Invalid annual bands now produce an informative variable/year-specific error.

Requires R runtime smoke testing and annual SWE regression checks; not executed here.
