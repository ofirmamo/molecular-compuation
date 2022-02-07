import random
from typing import List

from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from lab.common import *

random.seed(1)
fail_prob = 0.01


def create_primers(temp1: Dseqrecord, temp2: Dseqrecord):
    return [''.join(COMPLEMENTS[c] for c in temp1.seq.watson)[:10],
            ''.join(COMPLEMENTS[c] for c in temp2.seq.crick)[-10:]]


def pcr(tube: List[Dseqrecord], primers, rounds):
    for _ in range(rounds):
        tube = pcr_round(tube, primers)
    return tube


def anneal_and_amplify_watson(watson: str, primer: str):
    index = watson.find(''.join([COMPLEMENTS[c] for c in primer]))
    if index == 0:
        amplify = ''.join([COMPLEMENTS[c] for c in watson[len(primer):]])
        crick = primer + amplify
        return True, Dseqrecord(Dseq(watson=watson, crick=crick[::-1], ovhg=0, linear=True))
    return False, None


def anneal_and_amplify_crick(crick: str, primer: str):
    index = crick.rfind((''.join([COMPLEMENTS[c] for c in primer]))[::-1])
    if index == len(crick) - len(primer):
        amplify = ''.join([COMPLEMENTS[c] for c in crick[:-len(primer)]])
        watson = amplify + primer[::-1]
        return True, Dseqrecord(Dseq(watson=watson[::-1], crick=crick, ovhg=0, linear=True))
    return False, None


def heat_tube(tube: List[Dseqrecord]):
    result = []
    for molecule in tube:
        if molecule.seq.watson == '' or molecule.seq.crick == '' or random.random() < fail_prob:
            result.append(molecule)
            continue
        result.append(Dseqrecord(Dseq(watson=molecule.seq.watson, crick='', ovhg=0, linear=True)))
        result.append(Dseqrecord(Dseq(watson='', crick=molecule.seq.crick, ovhg=0, linear=True)))
    return result


def pcr_round(tube: [List[Dseqrecord]], primers: List[str]):
    result = []
    for molecule in tube:
        amplified = False
        for primer in primers:
            f, seq = (anneal_and_amplify_watson, molecule.seq.watson) if molecule.seq.watson != '' else (
                anneal_and_amplify_crick, molecule.seq.crick)
            success, amplified_molecule = f(seq, primer)
            if success and random.random() > fail_prob:
                amplified = True
                result.append(amplified_molecule)
                break
        if not amplified:
            result.append(molecule)
    return result
