"""Fill a driver configuration template with this host's paths.

    python render_driver_config.py TEMPLATE OUT key=value [key=value ...]

Only the named placeholders are replaced ({bin}, {out}, ...): the templates also carry the
drivers' own format fields, such as the beam width {w}, which must reach the driver as
written. Every value is made absolute, since the drivers resolve relative paths against
their own location. A placeholder left unfilled is an error.
"""

import json
import re
import sys
from pathlib import Path


def render(text, values):
    for key, value in values.items():
        text = text.replace("{" + key + "}", value)
    return text


def main():
    template, out, *pairs = sys.argv[1:]
    values = {}
    for pair in pairs:
        key, value = pair.split("=", 1)
        values[key] = str(Path(value).resolve())
    text = render(Path(template).read_text(), values)
    known = set(re.findall(r"\{([a-z0-9_]+)\}", Path(template).read_text()))
    # {w} and its kind are the drivers' own fields; anything else unfilled is a missing path.
    missing = sorted(k for k in known - set(values) if k not in {"w"})
    if missing:
        sys.exit(f"{template}: no value for {missing}")
    json.loads(text)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(text)


if __name__ == "__main__":
    main()
