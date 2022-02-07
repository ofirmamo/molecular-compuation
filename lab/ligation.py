import random
from typing import List

from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from lab.common import *


# Ligation is chemical reaction, therefore we will do it like a chemical reaction.
# First there is failure probability.
# Second, each round each fragment is picked up and chooses what to attach.

def ligate(tube: List[Dseqrecord], limit, rounds=7500):
    failure_probability = lambda molecule: len(molecule) / 10000
    random.seed(1)
    for i in range(rounds):
        molecule = random.choice(tube)
        p = failure_probability(molecule)
        if random.random() <= -1:
            print('Failure ligation, probability')
            continue
        else:
            ligation_single_fragment(molecule, tube, limit)


def end_watson_match(fragment: Dseqrecord, lookup_seq: str):
    diff = (len(fragment.seq.watson.strip()) + fragment.seq.ovhg) - len(fragment.seq.crick.strip())
    if diff < len(lookup_seq):
        return False
    return lookup_seq == fragment.seq.watson.strip()[-len(lookup_seq):]


def start_watson_match(fragment: Dseqrecord, lookup_seq: str):
    if -fragment.seq.ovhg < len(lookup_seq) and fragment.seq.crick.strip() != '':
        return False
    return lookup_seq == fragment.seq.watson.strip()[:len(lookup_seq)]


def end_crick_match(fragment: Dseqrecord, lookup_seq: str):
    diff = len(fragment.seq.crick.strip()) - (len(fragment.seq.watson.strip()) + fragment.seq.ovhg)
    if diff < len(lookup_seq):
        return False
    return lookup_seq == fragment.seq.crick.strip()[::-1][-len(lookup_seq):]


def start_crick_match(fragment: Dseqrecord, lookup_seq):
    if fragment.seq.ovhg < len(lookup_seq) and fragment.seq.watson.strip() != '':
        return False
    return lookup_seq == fragment.seq.crick.strip()[::-1][:len(lookup_seq)]


def ligation_single_watson_end(fragment: Dseqrecord, candidates: List[Dseqrecord], limit):
    lookup_seq = ''.join(
        [COMPLEMENTS[c] for c in fragment.seq.crick.strip()[::-1]]) if fragment.seq.watson.strip() == '' \
        else ''.join([COMPLEMENTS[c] for c in fragment.seq.crick.strip()[::-1][:fragment.seq.ovhg]])
    length = len(lookup_seq)
    for _ in range(length - limit + 1):
        for i, candidate in enumerate(candidates):
            if end_watson_match(candidate, lookup_seq):
                # del candidates[i]
                candidates.append(Dseqrecord(Dseq(
                    watson=(candidate.seq.watson.strip() + (
                            ' ' * (length - len(lookup_seq))) + fragment.seq.watson.strip()).strip(),
                    crick=(fragment.seq.crick.strip() + candidate.seq.crick.strip()).strip(),
                    ovhg=candidate.seq.ovhg if candidate.seq.crick.strip() != '' else -(
                            len(candidate.seq.watson.strip()) - len(lookup_seq)))))
                return True
        lookup_seq = lookup_seq[:-1]
    return False


def ligation_single_watson_start(fragment: Dseqrecord, candidates: List[Dseqrecord], limit):
    diff = len(fragment.seq.crick.strip()) - (len(fragment.seq.watson.strip()) + fragment.seq.ovhg)
    lookup_seq = ''.join([COMPLEMENTS[c] for c in fragment.seq.crick.strip()[::-1][-diff:]])
    length = len(lookup_seq)
    for _ in range(length - limit + 1):
        for i, candidate in enumerate(candidates):
            if start_watson_match(candidate, lookup_seq):
                # del candidates[i]
                candidates.append(Dseqrecord(Dseq(
                    watson=(fragment.seq.watson.strip() + (
                            ' ' * (length - len(lookup_seq))) + candidate.seq.watson.strip()).strip(),
                    crick=(candidate.seq.crick.strip() + fragment.seq.crick.strip()).strip(),
                    ovhg=fragment.seq.ovhg if fragment.seq.watson.strip() != '' else (
                            length - len(lookup_seq)))))
                return True
        lookup_seq = lookup_seq[1:]
    return False


def ligation_single_crick_end(fragment: Dseqrecord, candidates: List[Dseqrecord], limit):
    lookup_seq = ''.join([COMPLEMENTS[c] for c in fragment.seq.watson.strip()[:-fragment.seq.ovhg]])
    length = len(lookup_seq)
    for _ in range(length - limit + 1):
        for i, candidate in enumerate(candidates):
            if end_crick_match(candidate, lookup_seq):
                # del candidates[i]
                candidates.append(Dseqrecord(Dseq(
                    watson=(candidate.seq.watson.strip() + fragment.seq.watson.strip()).strip(),
                    crick=(fragment.seq.crick.strip() + (
                            ' ' * (length - len(lookup_seq))) + candidate.seq.crick.strip()).strip(),
                    ovhg=candidate.seq.ovhg if candidate.seq.watson != '' else len(candidate.seq.crick) - len(
                        lookup_seq))))
                return True
        lookup_seq = lookup_seq[:-1]
    return False


def ligation_single_crick_start(fragment: Dseqrecord, candidates: List[Dseqrecord], limit):
    diff = (len(fragment.seq.watson.strip()) + fragment.seq.ovhg) - len(fragment.seq.crick.strip())
    lookup_seq = ''.join([COMPLEMENTS[c] for c in fragment.seq.watson.strip()[-diff:]])
    length = len(lookup_seq)
    for _ in range(length - limit + 1):
        for i, candidate in enumerate(candidates):
            if start_crick_match(candidate, lookup_seq):
                # del candidates[i]
                candidates.append(Dseqrecord(Dseq(
                    watson=(fragment.seq.watson.strip() + candidate.seq.watson.strip()).strip(),
                    crick=(candidate.seq.crick.strip() + (
                            ' ' * (length - len(lookup_seq))) + fragment.seq.crick.strip()).strip(),
                    ovhg=fragment.seq.ovhg if candidate.seq.watson.strip() != '' and fragment.seq.crick.strip() != '' else - (
                            length - len(lookup_seq)))))
                return True
        lookup_seq = lookup_seq[1:]
    return False


def ligation_single_fragment(fragment: Dseqrecord, candidates: List[Dseqrecord], limit):
    case1_longer_crick_at_start = (
                                          fragment.seq.watson.strip() == '' and fragment.seq.crick.strip() != '') or fragment.seq.ovhg > 0
    case2_longer_watson_at_start = (
                                           fragment.seq.watson.strip() != '' and fragment.seq.crick.strip() == '') or fragment.seq.ovhg < 0
    case3_shorter_crick_at_end = len(fragment.seq.watson.strip()) + fragment.seq.ovhg > len(fragment.seq.crick.strip())
    case4_shorter_watson_at_end = len(fragment.seq.watson.strip()) + fragment.seq.ovhg < len(fragment.seq.crick.strip())

    if case4_shorter_watson_at_end and ligation_single_watson_start(fragment, candidates, limit):
        return True

    if case1_longer_crick_at_start and ligation_single_watson_end(fragment, candidates, limit):
        return True

    if case2_longer_watson_at_start and ligation_single_crick_end(fragment, candidates, limit):
        return True

    if case3_shorter_crick_at_end and ligation_single_crick_start(fragment, candidates, limit):
        return True

    return False


def should_ligate(fragment: Dseqrecord):
    return fragment.seq.ovhg != 0 or len(fragment.seq.watson.strip()) != len(fragment.seq.crick.strip()) or \
           (fragment.seq.watson.strip() != '' and fragment.seq.crick.strip() == '') or \
           (fragment.seq.watson.strip() == '' and fragment.seq.crick.strip() != '')
