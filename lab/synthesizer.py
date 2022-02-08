import random

from typing import List, Dict, Tuple
from pydna.dseq import Dseq as DnaSeq
from pydna.dseqrecord import Dseqrecord as DnaSeqRecord
from lab.common import base_list, COMPLEMENTS


def synthesize_watson_random(elements: List[str], bp_length) -> Dict[str, DnaSeqRecord]:
    res = dict()
    for element in elements:
        while True:
            dnaseq = DnaSeqRecord(
                DnaSeq(watson=''.join(random.choices(base_list, k=bp_length)), crick='', ovhg=0),
                linear=True,
                name=element)
            if valid_dna_strand(res, dnaseq, hamming_threshold=bp_length // 4):
                res[element] = dnaseq
                break
    return res


def join_watson(elements: List[Tuple[str, DnaSeqRecord, DnaSeqRecord]], take_first, take_second) -> Dict[str, DnaSeqRecord]:
    return {
        name: DnaSeqRecord(
            DnaSeq(watson=e1.seq.watson[-take_first:] + e2.seq.watson[:take_second], linear=True, crick='', ovhg=0),
            linear=True,
            name=name)
        for name, e1, e2 in elements
    }


def synthesize_crick_complements(elements) -> Dict[str, DnaSeqRecord]:
    return {
        element: DnaSeqRecord(
            DnaSeq(watson='', crick=''.join(COMPLEMENTS[c] for c in dna.seq.watson)[::-1], linear=True, ovhg=0))
        for element, dna in elements
    }


def valid_dna_strand(synthesized_elements_to_dna: Dict[str, DnaSeqRecord], dnaseq: DnaSeqRecord, hamming_threshold):
    for dna in synthesized_elements_to_dna.values():
        hamming_distance = sum(m1 != m2 for m1, m2 in zip(dna.seq.watson, dnaseq.seq.watson))
        if hamming_distance <= hamming_threshold:
            return False
    return True
