from typing import List

from pydna.dseqrecord import Dseqrecord


def magnetic_beads_filter(tube: List[Dseqrecord], seq: str):
    bead = lambda x: seq in x.seq.watson
    result = []
    for s in tube:
        if bead(s):
            result.append(s)
    return result
