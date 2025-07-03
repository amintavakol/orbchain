import random



from chemutils.Common.Util import (
    smi_to_unique_smi_fast,
    smi_to_unique_smi_map,
)

from atom.utils import *


def source_to_non_reactive_inter(ops, source, sink, samples):
    out = []
    for op in ops:
        if smi_to_unique_smi_map(op.srcAtom.connectedSmiles) != smi_to_unique_smi_map(
            source
        ):  # source must be there
            continue
        if smi_to_unique_smi_fast(
            op.sinkAtom.connectedSmiles
        ) == smi_to_unique_smi_fast(source):
            continue
        if smi_to_unique_smi_map(op.sinkAtom.connectedSmiles) == smi_to_unique_smi_map(
            sink
        ):  # don't want to have the actual sink
            continue

        out.append(op)

    if len(out) <= samples:
        return out
    else:
        return random.sample(out, samples)


def sink_to_non_reactive_inter(ops, source, sink, samples):
    out = []
    for op in ops:
        if smi_to_unique_smi_map(op.sinkAtom.connectedSmiles) != smi_to_unique_smi_map(
            sink
        ):  # sink must be there
            continue
        if smi_to_unique_smi_fast(op.srcAtom.connectedSmiles) == smi_to_unique_smi_fast(
            sink
        ):
            continue
        if smi_to_unique_smi_map(op.srcAtom.connectedSmiles) == smi_to_unique_smi_map(
            source
        ):  # don't want to have the actual source
            continue

        out.append(op)

    if len(out) <= samples:
        return out
    else:
        return random.sample(out, samples)


def non_reactive_pairs_inter(ops, source, sink, samples):
    out = []
    for op in ops:
        if smi_to_unique_smi_map(op.sinkAtom.connectedSmiles) == smi_to_unique_smi_map(
            sink
        ):  # should not be the actual sink
            continue
        if smi_to_unique_smi_map(op.srcAtom.connectedSmiles) == smi_to_unique_smi_map(
            source
        ):  # should not be the actual source
            continue
        if smi_to_unique_smi_fast(op.srcAtom.connectedSmiles) == smi_to_unique_smi_fast(
            op.srcAtom.connectedSmiles
        ):  # must be inter-molecular
            continue

        out.append(op)

    if len(out) <= samples:
        return out
    else:
        return random.sample(out, samples)


def non_reactive_pairs_intra(ops, source, sink, samples):
    out = []
    for op in ops:
        if smi_to_unique_smi_map(op.sinkAtom.connectedSmiles) == smi_to_unique_smi_map(
            sink
        ):  # should not be the actual sink
            continue
        if smi_to_unique_smi_map(op.srcAtom.connectedSmiles) == smi_to_unique_smi_map(
            source
        ):  # should not be the actual source
            continue
        if smi_to_unique_smi_fast(op.srcAtom.connectedSmiles) != smi_to_unique_smi_fast(
            op.srcAtom.connectedSmiles
        ):  # must be intra-molecular
            continue

        out.append(op)

    if len(out) <= samples:
        return out
    else:
        return random.sample(out, samples)


def source_to_non_reactive_intra(ops, source, sink, samples):
    out = []
    for op in ops:
        if smi_to_unique_smi_map(op.srcAtom.connectedSmiles) != smi_to_unique_smi_map(
            source
        ):  # source must be there
            continue
        if smi_to_unique_smi_fast(
            op.sinkAtom.connectedSmiles
        ) != smi_to_unique_smi_fast(source):
            continue
        if smi_to_unique_smi_map(op.sinkAtom.connectedSmiles) == smi_to_unique_smi_map(
            sink
        ):  # don't want to have the actual sink
            continue

        out.append(op)

    if len(out) <= samples:
        return out
    else:
        return random.sample(out, samples)


def sink_to_non_reactive_intra(ops, source, sink, samples):
    out = []
    for op in ops:
        if smi_to_unique_smi_map(op.sinkAtom.connectedSmiles) != smi_to_unique_smi_map(
            sink
        ):  # sink must be there
            continue
        if smi_to_unique_smi_fast(
            op.sinkAtom.connectedSmiles
        ) != smi_to_unique_smi_fast(source):
            continue
        if smi_to_unique_smi_map(op.srcAtom.connectedSmiles) == smi_to_unique_smi_map(
            source
        ):  # don't want to have the actual sink
            continue

        out.append(op)

    if len(out) <= samples:
        return out
    else:
        return random.sample(out, samples)
