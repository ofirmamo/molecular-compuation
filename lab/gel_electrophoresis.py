def gel_electrophoresis(tube, variant="all"):
    gel_container = sorted(tube, reverse=True, key=lambda x: (len(x.seq.watson), len(x.seq.crick)))
    if variant == "max":
        return gel_container[0] if any(gel_container) else None
    if variant == "min":
        return gel_container[-1] if any(gel_container) else None
    return list(gel_container)  # variant == "all"
