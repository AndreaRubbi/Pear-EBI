#


### parser
[source](https://github.com/AndreaRubbi/Pear-EBI/blob/pear_ebi/pear_ebi/tree_emb_parser.py/#L6)
```python
.parser(
   argv = None
)
```

---
Parses PEAR's command line arguments.


**Args**

* **argv** (list, optional) : argument list to parse. Defaults to None, which means
    sys.argv[1:], i.e. the real command line. Passing a list makes the parser
    testable and lets callers drive PEAR programmatically; without it there
    was no way to exercise this function at all.


**Returns**

* **Namespace**  : the parsed arguments
