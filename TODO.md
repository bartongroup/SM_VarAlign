# TODO

- Refactor or remove unnessential results files, e.g. UMD/CME residue lists, etc.
- Put temp test data in better location, annoying when note cleared
- Better results file organisation, e.g. no hidden folders
- Review exceptions
- Remove odd patterns, e.g. "action" argument in interpret regression
- Remove non-essential stuff, e.g. interpret regression
- Review and streamline logging
- Review all output files
- Consider using mypy pre-commit hooks


## Notes

### mypy pre-commit conf

```python
-   repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.11.0
    hooks:
    - id: mypy
      args: ["--install-types", "--non-interactive"]
```
