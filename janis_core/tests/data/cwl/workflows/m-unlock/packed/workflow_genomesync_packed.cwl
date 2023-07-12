{
    "class": "CommandLineTool",
    "label": "Genome synchronisation with EBI",
    "doc": "Runs the genome synchronisation application which downloads genomes from the EBI into the corresponding taxonomic path\n",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [
        {
            "type": [
                "null",
                "int"
            ],
            "doc": "Maximum number of genomes per bin",
            "label": "bin size",
            "inputBinding": {
                "prefix": "-bin"
            },
            "default": 10,
            "id": "#main/bin"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Filtering by taxonomic class name",
            "label": "class name",
            "inputBinding": {
                "prefix": "-class"
            },
            "id": "#main/clazz"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "label": "Output Destination",
            "doc": "Optional Output destination used for cwl-prov reporting.",
            "id": "#main/destination"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "doc": "Retrieves genome files in embl format",
            "label": "embl format",
            "inputBinding": {
                "prefix": "-embl"
            },
            "id": "#main/embl"
        },
        {
            "type": "string",
            "doc": "Path of the local enaBrowserTools git repository",
            "label": "enaBrowserTools location",
            "inputBinding": {
                "prefix": "-enaBrowserTools"
            },
            "default": "/unlock/infrastructure/binaries/enaBrowserTools",
            "id": "#main/enaBrowserTools"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Filtering by taxonomic family name",
            "label": "family name",
            "inputBinding": {
                "prefix": "-family"
            },
            "id": "#main/family"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "doc": "Retrieves genome files in fasta format",
            "label": "fasta format",
            "inputBinding": {
                "prefix": "-fasta"
            },
            "id": "#main/fasta"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Filtering by taxonomic genus name",
            "label": "genus name",
            "inputBinding": {
                "prefix": "-genus"
            },
            "id": "#main/genus"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "doc": "ensures obtained data is stored in irods",
            "label": "irods connection",
            "inputBinding": {
                "prefix": "-irods"
            },
            "id": "#main/irods"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "doc": "Keep the downloaded files from ENA after successful completion",
            "label": "downloaded files are stored permanently",
            "inputBinding": {
                "prefix": "-keep"
            },
            "id": "#main/keep"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Filtering by taxonomic order name",
            "label": "order name",
            "inputBinding": {
                "prefix": "-order"
            },
            "id": "#main/order"
        },
        {
            "type": "string",
            "doc": "output folder to store the results in",
            "label": "output folder",
            "inputBinding": {
                "prefix": "-output"
            },
            "id": "#main/output"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Filtering by taxonomic phylum name",
            "label": "phylum name",
            "inputBinding": {
                "prefix": "-phylum"
            },
            "id": "#main/phylum"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Skips the conversion of the strains to rdf",
            "label": "skip rdf conversion",
            "inputBinding": {
                "prefix": "-skip"
            },
            "id": "#main/skip"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Filtering by taxonomic species name",
            "label": "species name",
            "inputBinding": {
                "prefix": "-species"
            },
            "id": "#main/species"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Enable strict mode, stops completelty when a conversion fails",
            "label": "strict mode",
            "inputBinding": {
                "prefix": "-strict"
            },
            "id": "#main/strict"
        },
        {
            "type": "string",
            "doc": "input of the assembly_summary_genbank file",
            "label": "assembly_summary_genbank file",
            "inputBinding": {
                "prefix": "-summary"
            },
            "default": "assembly_summary_genbank.txt",
            "id": "#main/summary"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "doc": "Filtering by taxonomic superkingdom name",
            "label": "superkingdom name",
            "inputBinding": {
                "prefix": "-superkingdom"
            },
            "id": "#main/superkingdom"
        },
        {
            "type": "string",
            "doc": "Path of the taxonomy.hdt lookup file",
            "label": "Taxnomy RDF (HDT) file",
            "inputBinding": {
                "prefix": "-taxonomy"
            },
            "default": "/tempZone/References/Databases/UniProt/taxonomy.hdt",
            "id": "#main/taxonomy"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "doc": "Force update of the assembly_summary_genbank lookup list",
            "label": "force lookup list update",
            "inputBinding": {
                "prefix": "-update"
            },
            "id": "#main/update"
        }
    ],
    "arguments": [
        "java",
        "-Xmx3g",
        "-jar",
        "/unlock/infrastructure/binaries/GenomeSync.jar"
    ],
    "id": "#main",
    "https://schema.org/author": [
        {
            "class": "https://schema.org/Person",
            "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
            "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
            "https://schema.org/name": "Jasper Koehorst"
        }
    ],
    "https://schema.org/citation": "https://m-unlock.nl",
    "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
    "https://schema.org/dateCreated": "2020-00-00",
    "https://schema.org/dateModified": "2022-05-00",
    "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
    "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential",
    "outputs": [
        {
            "type": "Directory",
            "outputBinding": {
                "glob": "bla"
            },
            "id": "#main/test"
        }
    ],
    "cwlVersion": "v1.2",
    "$namespaces": {
        "s": "https://schema.org/"
    }
}
