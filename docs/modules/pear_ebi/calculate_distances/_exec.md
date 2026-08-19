#


## PearExecutableError
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/_exec.py/#L41)
```python 
PearExecutableError()
```


---
A bundled native tool could not be run, or failed while running.

----


### run_process
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/_exec.py/#L97)
```python
.run_process(
   cmd_list, *, timeout = None
)
```

---
Run a native tool, always capturing both streams.

Never raises for a non-zero exit; inspect the returned CompletedRun. Failures
that happen before exec are mapped onto the RC_* sentinels so callers can give
a specific message.

----


### resolve_binary
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/_exec.py/#L143)
```python
.resolve_binary(
   system_name, packaged_path, *, tool_label
)
```

---
Pick the executable to use, preferring one already on PATH.

A binary on PATH wins over the bundled copy, which is how a locally built tool
is honoured -- but it used to happen silently, so a stale or incompatible
``hashrf`` on PATH would be used with no indication. That now warns.

Raises PearExecutableError if neither is available.

----


### remove_file
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/calculate_distances/_exec.py/#L242)
```python
.remove_file(
   path
)
```

---
Delete a file, reporting failure instead of hiding it.

Replaces ``bash_command(f"rm {path}")``, which shelled out with an unquoted
path (so it broke on spaces), sent both streams to DEVNULL, and returned 0
unconditionally whether or not the removal worked.
