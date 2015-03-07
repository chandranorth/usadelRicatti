def mangle_docstrings(app, what, name, obj, options, lines):
    new_lines = []
    skip_empty = False
    indent = ""
    for j, line in enumerate(lines):
        if skip_empty and not line.strip():
            continue
        else:
            skip_empty = False

        try:
            if lines[j+1].startswith('----'):
                if line.strip() in ('Parameters','Returns','Yields','Raises',
                                    'Attributes'):
                    new_lines.append(":%s:" % line)
                    indent = "    "
                else:
                    new_lines.append(".. rubric:: %s" % line)
                    new_lines.append("")
                    indent = ""
                skip_empty = True
                lines[j+1] = ""
                continue
        except IndexError:
            pass
        new_lines.append(indent + line)

    lines[:] = new_lines

def setup(app):
    app.connect('autodoc-process-docstring', mangle_docstrings)

