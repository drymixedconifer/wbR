# All-NA output band handling

`wb_atomic_raster()` now permits entirely missing bands only when they were also entirely missing in its input. Band counts, raster geometry and readable cell counts are checked before atomic replacement. All-NA bands trigger warnings with band names, not zeros. A write that loses non-NA source values still fails.

This allows first-year seasonal snow bands to be represented as missing without misreporting zero snow. This does not resolve missing seasonal history in the scientific output.
