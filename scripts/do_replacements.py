
from typing import Tuple
import argparse
import regex as re
from dataclasses import dataclass
from copy import deepcopy


def main():
    args = parse_args()
    template_f, replacements_f = read_files(args)
    template_segments = gather_segments(template_f)
    replacement_segments = gather_segments(replacements_f)
    outfile = replace_segments(template_f, template_segments, replacement_segments)
    write_file(args, outfile)

### GATHERING SEGMENTS

@dataclass
class Segment:
    label: str
    start: int
    end: int
    contents: str

segment_matcher = r'\[\/\/\]: <> \(([^/].*?)\)([\s\S]*?)\[\/\/\]: <> \(([/].*?)\)'

def gather_segments(text: str) -> dict[str, Segment]:
    segments: dict[str, Segment] = {}
    iterator = re.finditer(segment_matcher, text)
    matches = [match for match in iterator]
    for m in matches:
        s = Segment(
            label=m.group(1),
            start=m.start(),
            end=m.end(),
            contents=m.group(2)
        )
        segments[s.label] = s
    return segments

### REPLACING SEGMENTS

def replace_segments(template: str, template_segments: dict[str, Segment], replacement_segments: dict[str, Segment]) -> str:
    output: str = deepcopy(template)
    tsegs = sorted(template_segments.values(), key=lambda s: s.start, reverse=True)
    for tseg in tsegs:
        rseg = replacement_segments[tseg.label]
        output = output[:tseg.start] + rseg.contents + output[tseg.end:]
    return output

### FILE IO

def read_files(args: argparse.Namespace) -> Tuple[str, str]:
    template = open(args.infile, 'r').read()
    replacements = open(args.replacements, 'r').read()
    return template, replacements

def write_file(args: argparse.Namespace, outfile: str):
    with open(args.output, 'w') as f:
        f.write(outfile)

### ARGS

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str, help='input filepath')
    parser.add_argument('-r', '--replacements', type=str, help='file containing replacements')
    parser.add_argument('-o', '--output', type=str, default='replaced.md', help='output filepath') 
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()