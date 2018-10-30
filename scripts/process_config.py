#!/usr/bin/env python
import yaml
import sys


def main():

    base_config = sys.argv[1]
    out_config = sys.argv[2]
    config = sys.argv[3:]

    # sys.argv ignores quotes, so would split e.g.
    # MINI_ASSEMBLE_OPTS="-n 10" GUPPY_OPTS="--hp_correct 1"
    # into
    # ['MINI_ASSEMBLE_OPTS=-n', '10', 'GUPPY_OPTS=--hp_correct', '1']

    fixed = []
    for c in config:
        if '=' in c:
            fixed.append(c)
        else:
            fixed[-1] += ' ' + c
    config = fixed

    conf = yaml.load(open(base_config))

    d = {}
    for arg in config:
        k, v = arg.split('=')
        d[k] = v

    conf.update(d)

    yaml.dump(conf, open(out_config, 'w'))


if __name__ == '__main__':
    main()
