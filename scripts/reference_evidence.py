"""Evidence vocabulary shared by taxon priors and repeat interpretation.

Membership permits a supported reference annotation; it does not establish
terminal localization in the current sample. Reference status and taxonomic
scope remain separate gates, including candidate and unknown exceptions.
"""

SUPPORTED_REFERENCE_EVIDENCE = frozenset({
    'experimental', 'experimental_terminal', 'assembly_terminal',
    'published_terminal_long_read', 'long_read_terminal',
    'published_species_reference', 'clade_consensus',
})
