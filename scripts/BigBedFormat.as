table ABP
"collapses from ABP"
    (
    string chrom;               "Reference sequence chromosome or scaffold"
    uint   chromStart;          "Start position in chromosome"
    uint   chromEnd;            "End position in chromosome"
    string name;                "Name or ID of item, ideally both human-readable and unique"
    uint score;                 "Score (0-1000)"
    char[1] strand;             "+ or - for strand"
    uint thickStart;            "Start of where display should be thick (start codon)"
    uint thickEnd;              "End of where display should be thick (stop codon)"
    uint reserved;              "RGB value (use R,G,B string in input file)"
    string collapse;            "Collapse name"
        string description;             "description of the collapse"
    )
