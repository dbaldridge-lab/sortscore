"""
Variant-string parsing helpers.
"""
import re


def is_aa_change_format(variant_seq: str) -> bool:
    """
    Detect whether a sequence uses pre-annotated AA-change notation.
    """
    clean_seq = re.sub(r'^p\.', '', variant_seq.strip())
    clean_seq = re.sub(r'[\s\.\-_]', '', clean_seq)

    aa_codes = [
        'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
        'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',
        'Ter', 'Amb', 'Sec', 'Pyl'
    ]

    aa_pattern = '|'.join(aa_codes)
    pattern_3letter = rf'^({aa_pattern})([0-9]+)({aa_pattern})$'
    if re.match(pattern_3letter, clean_seq, re.IGNORECASE):
        return True

    pattern_1letter = r'^[A-Z*X][0-9]+[A-Z*X]$'
    if re.match(pattern_1letter, clean_seq, re.IGNORECASE):
        return True

    return False


def parse_aa_change(variant_seq: str) -> tuple:
    """
    Parse a pre-annotated AA change into (ref_aa, position, alt_aa).
    """
    clean_seq = re.sub(r'^p\.', '', variant_seq.strip())
    clean_seq = re.sub(r'[\s\.\-_]', '', clean_seq)

    aa_map = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
        'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
        'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
        'Tyr': 'Y', 'Val': 'V', 'Ter': '*', 'Amb': 'X', 'Sec': 'U', 'Pyl': 'O'
    }

    aa_codes = list(aa_map.keys())
    aa_pattern = '|'.join(aa_codes)

    pattern = rf'^({aa_pattern})([0-9]+)({aa_pattern})$'
    match = re.match(pattern, clean_seq, re.IGNORECASE)
    if match:
        ref_aa_3 = match.group(1).capitalize()
        position = int(match.group(2))
        alt_aa_3 = match.group(3).capitalize()

        ref_aa = aa_map.get(ref_aa_3, ref_aa_3[0])
        alt_aa = aa_map.get(alt_aa_3, alt_aa_3[0])
        return ref_aa, position, alt_aa

    pattern = r'^([A-Z*X])([0-9]+)([A-Z*X])$'
    match = re.match(pattern, clean_seq, re.IGNORECASE)
    if match:
        ref_aa = match.group(1).upper()
        position = int(match.group(2))
        alt_aa = match.group(3).upper()
        return ref_aa, position, alt_aa

    raise ValueError(f"Invalid AA change format: {variant_seq}")
