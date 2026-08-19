#


### platform_bin_dir
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/_install_helpers.py/#L45)
```python
.platform_bin_dir(
   package_root = None
)
```

---
Return the per-platform binary root for this OS, or None if unsupported.

Returns e.g. .../calculate_distances/linux_bin. The directory is returned even
when it does not exist, so callers can produce a useful message; use
``os.path.isdir`` to test. Returns None on platforms we ship no binaries for
(notably Windows).

----


### describe_platform
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/_install_helpers.py/#L59)
```python
.describe_platform()
```

---
Short human-readable platform string, used in error messages.

Includes the machine architecture because the bundled macOS binaries are
arm64-only: an Intel Mac fails with an architecture error, not a missing file.

----


### hashrf_binary
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/_install_helpers.py/#L68)
```python
.hashrf_binary(
   package_root = None
)
```

---
Path to the bundled hashrf executable for this platform, or None.

----


### tqdist_bin_dir
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/_install_helpers.py/#L76)
```python
.tqdist_bin_dir(
   package_root = None
)
```

---
Directory holding the bundled tqDist executables, or None.

----


### native_executables
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/_install_helpers.py/#L84)
```python
.native_executables(
   package_root = None
)
```

---
Every bundled native executable that should carry the exec bit.

----


### ensure_native_executables
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/_install_helpers.py/#L99)
```python
.ensure_native_executables(
   package_root = None
)
```

---
Ensure the bundled native executables are executable on POSIX systems.

Safe to call at import time: no-ops on non-POSIX platforms and on platforms we
ship no binaries for, and skips files that are already executable or absent.

Returns the list of paths whose mode was changed. Raises nothing; problems are
reported through ``warnings.warn`` so a read-only or shared install is visible
without breaking the import.

----


### build_tqdist
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/_install_helpers.py/#L142)
```python
.build_tqdist(
   package_root = None, use_cmake_first = True, timeout = 300
)
```

---
Attempt to (re)build the bundled tqDist native tools in place.

Intended for the case where a shipped binary will not run on the current
platform, for example because it was built for a different architecture
(the bundled macOS binaries are arm64-only).

Returns (success: bool, message: str).
