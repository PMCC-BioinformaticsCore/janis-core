{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "compress a file multithreaded with pigz",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/pigz:2.6",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.6"
                            ],
                            "specs": [
                                "https://anaconda.org/conda-forge/pigz"
                            ],
                            "package": "pigz"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "pigz",
                "-c"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.inputfile)"
                }
            ],
            "stdout": "$(inputs.inputfile.basename).gz",
            "inputs": [
                {
                    "type": "File",
                    "id": "#pigz.cwl/inputfile"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 1,
                    "inputBinding": {
                        "prefix": "-p"
                    },
                    "id": "#pigz.cwl/threads"
                }
            ],
            "id": "#pigz.cwl",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputfile.basename).gz"
                    },
                    "id": "#pigz.cwl/outfile"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "bash",
                "-x",
                "script.sh"
            ],
            "label": "CarveMe GEMstats",
            "doc": "Small summary of a list of CarveMe genome-scale metabolic models in sbml-fbc2 format\nContains; number of mets,reactions and genes\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "gemstats_output",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nidentifier=$1\nshift;\necho \"Model Mets Reactions Genes\" > $identifier\\_CarveMe_GEMstats.tsv\nfor file in \"$@\"\ndo\n  bash /unlock/infrastructure/scripts/GEMstats.sh $file\ndone >> $identifier\\_CarveMe_GEMstats.tsv"
                        }
                    ]
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "label": "CarveMe GEMs",
                    "doc": "List of CarveMe metabolic models in sbml-fbc2 format.",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#GEMstats.cwl/carveme_gems"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#GEMstats.cwl/identifier"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_CarveMe_GEMstats.tsv"
                    },
                    "id": "#GEMstats.cwl/carveme_GEMstats"
                }
            ],
            "id": "#GEMstats.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2022-01-01",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "carve"
            ],
            "label": "CarveMe",
            "doc": "CarveMe is a python-based tool for genome-scale metabolic model reconstruction.\n(Workflow will quit as successful even though no model can be created. Check messages.)\n    \n",
            "requirements": [
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "carve_output",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nsource /root/miniconda/bin/activate\nconda init bash\nconda activate /unlock/infrastructure/conda/carveme/cplex/carveme_1.5.1\ncarve $@"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "packages": [
                        {
                            "version": [
                                "1.5.1"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/carveme"
                            ],
                            "package": "carveme"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "outputs": [
                {
                    "label": "CarveMe GEM",
                    "doc": "CarveMe metabolic model Output SBML in sbml-fbc2 format",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.protein_file.nameroot).GEM.xml"
                    },
                    "id": "#carveme.cwl/carveme_gem"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Gap fill",
                    "doc": "Gap fill model for given media",
                    "inputBinding": {
                        "prefix": "--gapfill"
                    },
                    "id": "#carveme.cwl/gapfill"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Initial media",
                    "doc": "Initialize model with given medium",
                    "inputBinding": {
                        "prefix": "--init"
                    },
                    "id": "#carveme.cwl/init"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "inputBinding": {
                        "prefix": "--mediadb"
                    },
                    "id": "#carveme.cwl/mediadb"
                },
                {
                    "type": "File",
                    "label": "Input fasta file",
                    "doc": "Proteins sequence file in FASTA format.",
                    "inputBinding": {
                        "position": 0
                    },
                    "id": "#carveme.cwl/protein_file"
                }
            ],
            "arguments": [
                "--fbc2",
                {
                    "prefix": "--output",
                    "valueFrom": "$(inputs.protein_file.nameroot).GEM.xml"
                }
            ],
            "id": "#carveme.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2022-06-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "ExpressionTool",
            "doc": "Transforms the input files to a mentioned directory\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "id": "#files_to_folder.cwl/destination"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#files_to_folder.cwl/files"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "Directory"
                        }
                    ],
                    "id": "#files_to_folder.cwl/folders"
                }
            ],
            "expression": "${\n  var array = []\n  if (inputs.files != null) {\n    array = array.concat(inputs.files)\n  }\n  if (inputs.folders != null) {\n    array = array.concat(inputs.folders)\n  }\n  var r = {\n     'results':\n       { \"class\": \"Directory\",\n         \"basename\": inputs.destination,\n         \"listing\": array\n       } \n     };\n   return r; \n }\n",
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#files_to_folder.cwl/results"
                }
            ],
            "id": "#files_to_folder.cwl",
            "http://schema.org/citation": "https://m-unlock.nl",
            "http://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "http://schema.org/dateCreated": "2020-00-00",
            "http://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "http://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "MEMOTE",
            "doc": "MEMOTE, short for metabolic model testing",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "MEMOTE",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nsource /root/miniconda/bin/activate\nconda init bash\nconda activate /unlock/infrastructure/conda/memote/cplex/memote_0.13.0\nmemote $@"
                        }
                    ]
                }
            ],
            "baseCommand": [
                "bash",
                "script.sh"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Metabolic model (xml format)",
                    "label": "Metabolic model",
                    "inputBinding": {
                        "position": 100
                    },
                    "id": "#memote.cwl/GEM"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Report",
                    "doc": "Take a snapshot of a model's state and generate a report.",
                    "inputBinding": {
                        "position": 0
                    },
                    "id": "#memote.cwl/report_snapshot"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Run",
                    "doc": "Run the test suite on a single model and collect results.",
                    "inputBinding": {
                        "position": 0
                    },
                    "id": "#memote.cwl/run"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; find incorrect thermodynamic reversibility",
                    "label": "Skip test; find incorrect thermodynamic reversibility",
                    "inputBinding": {
                        "prefix": "--skip skip_test_find_incorrect_thermodynamic_reversibility",
                        "position": 15
                    },
                    "id": "#memote.cwl/skip_test_find_incorrect_thermodynamic_reversibility"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; Find metabolites consumed with closed bounds",
                    "label": "Skip test; Find metabolites consumed with closed bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_consumed_with_closed_bounds",
                        "position": 12
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_consumed_with_closed_bounds"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; Find metabolites not consumedwith open_bounds",
                    "label": "Skip test; Find metabolites not consumedwith open_bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_not_consumed_with_open_bounds",
                        "position": 14
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_not_consumed_with_open_bounds"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; Find metabolites not produced with open bounds",
                    "label": "Skip test; Find metabolites not produced with open bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_not_produced_with_open_bounds",
                        "position": 13
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_not_produced_with_open_bounds"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; find metabolites produced with closed bounds",
                    "label": "Skip test; find metabolites produced with closed bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_produced_with_closed_bounds",
                        "position": 11
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_produced_with_closed_bounds"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "solver",
                    "doc": "Set the solver to be used. [cplex|glpk|gurobi|glpk_exact]. default; glpk",
                    "inputBinding": {
                        "prefix": "--solver",
                        "position": 3
                    },
                    "id": "#memote.cwl/solver"
                }
            ],
            "arguments": [
                "${\n  if (inputs.run){\n    return \"run --filename \" + inputs.GEM.basename + \"_MEMOTE.json.gz\";\n  } else {\n    return '';\n  }\n}\n",
                "${\n  if (inputs.report_snapshot){\n    return \"report snapshot --filename \" + inputs.GEM.basename + \"_MEMOTE.html\";\n  } else {\n    return '';\n  }\n}\n"
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.GEM.basename)_MEMOTE.html"
                    },
                    "id": "#memote.cwl/report_html"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.GEM.basename)_MEMOTE.json.gz"
                    },
                    "id": "#memote.cwl/run_json"
                }
            ],
            "id": "#memote.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2022-06-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "Prodigal",
            "doc": "Prokaryotic gene prediction using Prodigal",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/prodigal:2.6.3",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.6.3"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/prodigal"
                            ],
                            "package": "prodigal"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "prodigal"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.input_fasta.nameroot).prodigal",
                    "prefix": "-o"
                },
                {
                    "valueFrom": "$(inputs.input_fasta.nameroot).prodigal.ffn",
                    "prefix": "-d"
                },
                {
                    "valueFrom": "$(inputs.input_fasta.nameroot).prodigal.faa",
                    "prefix": "-a"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.input_fasta.nameroot).prodigal.faa"
                    },
                    "id": "#prodigal.cwl/predicted_proteins_faa"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.input_fasta.nameroot).prodigal.ffn"
                    },
                    "id": "#prodigal.cwl/predicted_proteins_ffn"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.input_fasta.nameroot).prodigal"
                    },
                    "id": "#prodigal.cwl/predicted_proteins_out"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "-i"
                    },
                    "id": "#prodigal.cwl/input_fasta"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Input is a meta-genome",
                    "inputBinding": {
                        "prefix": "-p",
                        "valueFrom": "meta"
                    },
                    "id": "#prodigal.cwl/meta_mode"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Input is an isolate genome",
                    "inputBinding": {
                        "prefix": "-p",
                        "valueFrom": "single"
                    },
                    "id": "#prodigal.cwl/single_mode"
                }
            ],
            "id": "#prodigal.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-6867-2039",
                    "https://schema.org/name": "Ekaterina Sakharova"
                }
            ],
            "https://schema.org/copyrightHolder'": "EMBL - European Bioinformatics Institute",
            "https://schema.org/license'": "https://www.apache.org/licenses/LICENSE-2.0",
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2022-06-00",
            "https://schema.org/dateModified": "2022-08-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential",
            "https://schema.org/copyrightNotice": " Copyright < 2022 EMBL - European Bioinformatics Institute This file has been modified by UNLOCK - Unlocking Microbial Potential "
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "bash",
                "script.sh"
            ],
            "label": "SMETANA",
            "doc": "Species METabolic interaction ANAlysis (SMETANA) is a python-based command line tool to analyse microbial communities.\nIt takes as input a microbial communtity (from a collection of genome-scale metabolic models in SBML format) and \ncomputes of several metrics that describe the potential for cross-feeding interactions between community members.\n    \n",
            "requirements": [
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "smetana_output",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nsource /root/miniconda/bin/activate\nconda init bash\nconda activate /unlock/infrastructure/conda/smetana/cplex/smetana_1.2.0\nsmetana $@"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "label": "SMETANA output",
                    "doc": "SMETANA detailed tsv output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SMETANA*"
                    },
                    "id": "#smetana.cwl/detailed_output_tsv"
                }
            ],
            "inputs": [
                {
                    "label": "Metabolic model",
                    "doc": "Multiple Metabolic models (xml format)",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 10
                    },
                    "id": "#smetana.cwl/GEM"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Runs MIP/MRO and is much faster, recommended when analysing multiple communities.",
                    "label": "Detailed mode",
                    "inputBinding": {
                        "prefix": "--detailed",
                        "position": 11
                    },
                    "id": "#smetana.cwl/detailed"
                },
                {
                    "type": "string",
                    "label": "Flavor",
                    "doc": "Expected SBML flavor of the input files (cobra or fbc2)",
                    "inputBinding": {
                        "prefix": "--flavor",
                        "position": 3
                    },
                    "default": "fbc2",
                    "id": "#smetana.cwl/flavor"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Calculates all inter-species interactions (much slower), check Algorithms for details.",
                    "label": "Global mode",
                    "inputBinding": {
                        "prefix": "--global",
                        "position": 11
                    },
                    "default": false,
                    "id": "#smetana.cwl/global"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "Identifier used",
                    "id": "#smetana.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "inputBinding": {
                        "prefix": "-m",
                        "position": 1
                    },
                    "id": "#smetana.cwl/media"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "inputBinding": {
                        "prefix": "--mediadb",
                        "position": 2
                    },
                    "id": "#smetana.cwl/mediadb"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "solver",
                    "doc": "Set the solver to be used. [cplex|glpk|gurobi|glpk_exact]. default; glpk",
                    "inputBinding": {
                        "prefix": "--solver",
                        "position": 4
                    },
                    "id": "#smetana.cwl/solver"
                }
            ],
            "arguments": [
                "--verbose",
                {
                    "prefix": "--output",
                    "valueFrom": "$(inputs.identifier)_SMETANA"
                }
            ],
            "id": "#smetana.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2022-06-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "label": "Metagenomic GEM construction from assembly",
            "doc": "Workflow for Metagenomics from bins to metabolic model.<br>\nSummary\n  - Prodigal gene prediction\n  - CarveMe genome scale metabolic model reconstruction\n  - MEMOTE for metabolic model testing\n  - SMETANA Species METabolic interaction ANAlysis\n\nOther UNLOCK workflows on WorkflowHub: https://workflowhub.eu/projects/16/workflows?view=default<br><br>\n\n**All tool CWL files and other workflows can be found here:**<br>\n  Tools: https://gitlab.com/m-unlock/cwl<br>\n  Workflows: https://gitlab.com/m-unlock/cwl/workflows<br>\n\n**How to setup and use an UNLOCK workflow:**<br>\nhttps://m-unlock.gitlab.io/docs/setup/setup.html<br>\n",
            "outputs": [
                {
                    "label": "CarveMe GEMs folder",
                    "doc": "CarveMe metabolic models folder",
                    "type": "Directory",
                    "outputSource": "#main/carveme_files_to_folder/results",
                    "id": "#main/carveme_gems_folder"
                },
                {
                    "label": "GEMstats",
                    "doc": "CarveMe GEM statistics",
                    "type": "File",
                    "outputSource": "#main/gemstats/carveme_GEMstats",
                    "id": "#main/gemstats_out"
                },
                {
                    "label": "MEMOTE outputs folder",
                    "doc": "MEMOTE outputs folder",
                    "type": "Directory",
                    "outputSource": "#main/memote_files_to_folder/results",
                    "id": "#main/memote_folder"
                },
                {
                    "label": "Protein files folder",
                    "doc": "Prodigal predicted proteins (compressed) fasta files",
                    "type": "Directory",
                    "outputSource": "#main/prodigal_files_to_folder/results",
                    "id": "#main/protein_fasta_folder"
                },
                {
                    "label": "SMETANA output",
                    "doc": "SMETANA detailed output table",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/smetana/detailed_output_tsv",
                    "id": "#main/smetana_output"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Bin/genome fasta files",
                    "label": "Genome/bin",
                    "id": "#main/bins"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output Destination (prov only)",
                    "doc": "Not used in this workflow. Output destination used for cwl-prov reporting only.",
                    "id": "#main/destination"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Gap fill",
                    "doc": "Gap fill model for given media",
                    "id": "#main/gapfill"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "Identifier used",
                    "id": "#main/identifier"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "id": "#main/mediadb"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Run SMETANA",
                    "doc": "Run SMETANA (Species METabolic interaction ANAlysis)",
                    "default": false,
                    "id": "#main/run_smetana"
                },
                {
                    "type": "string",
                    "doc": "Solver to be used in MEMOTE and SMETANA (defaul; cplex)",
                    "default": "cplex",
                    "id": "#main/solver"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "number of threads",
                    "default": 2,
                    "id": "#main/threads"
                }
            ],
            "steps": [
                {
                    "label": "CarveMe",
                    "doc": "Genome-scale metabolic models reconstruction with CarveMe",
                    "run": "#carveme.cwl",
                    "scatter": [
                        "#main/carveme/protein_file"
                    ],
                    "in": [
                        {
                            "source": "#main/gapfill",
                            "id": "#main/carveme/gapfill"
                        },
                        {
                            "source": "#main/mediadb",
                            "id": "#main/carveme/mediadb"
                        },
                        {
                            "source": "#main/prodigal/predicted_proteins_faa",
                            "id": "#main/carveme/protein_file"
                        }
                    ],
                    "out": [
                        "#main/carveme/carveme_gem"
                    ],
                    "id": "#main/carveme"
                },
                {
                    "doc": "Preparation of workflow output files to a specific output folder",
                    "label": "CarveMe GEMs to folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "CarveMe_GEMs",
                            "id": "#main/carveme_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/compress_carveme/outfile"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/carveme_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/carveme_files_to_folder/results"
                    ],
                    "id": "#main/carveme_files_to_folder"
                },
                {
                    "label": "Compress GEM",
                    "doc": "Compress CarveMe GEM",
                    "run": "#pigz.cwl",
                    "scatter": "#main/compress_carveme/inputfile",
                    "in": [
                        {
                            "source": [
                                "#main/carveme/carveme_gem"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/compress_carveme/inputfile"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/compress_carveme/threads"
                        }
                    ],
                    "out": [
                        "#main/compress_carveme/outfile"
                    ],
                    "id": "#main/compress_carveme"
                },
                {
                    "label": "Compress proteins",
                    "doc": "Compress prodigal protein files",
                    "run": "#pigz.cwl",
                    "scatter": "#main/compress_prodigal/inputfile",
                    "in": [
                        {
                            "source": [
                                "#main/prodigal/predicted_proteins_faa"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/compress_prodigal/inputfile"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/compress_prodigal/threads"
                        }
                    ],
                    "out": [
                        "#main/compress_prodigal/outfile"
                    ],
                    "id": "#main/compress_prodigal"
                },
                {
                    "label": "GEM stats",
                    "doc": "CarveMe GEM statistics",
                    "run": "#GEMstats.cwl",
                    "in": [
                        {
                            "source": "#main/carveme/carveme_gem",
                            "id": "#main/gemstats/carveme_gems"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/gemstats/identifier"
                        }
                    ],
                    "out": [
                        "#main/gemstats/carveme_GEMstats"
                    ],
                    "id": "#main/gemstats"
                },
                {
                    "doc": "Preparation of workflow output files to a specific output folder",
                    "label": "MEMOTE output",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "MEMOTE",
                            "id": "#main/memote_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/memote_report_snapshot/report_html",
                                "#main/memote_run/run_json"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/memote_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/memote_files_to_folder/results"
                    ],
                    "id": "#main/memote_files_to_folder"
                },
                {
                    "label": "MEMOTE report snapshot",
                    "doc": "Take a snapshot of a model's state and generate a report.",
                    "run": "#memote.cwl",
                    "scatter": [
                        "#main/memote_report_snapshot/GEM"
                    ],
                    "in": [
                        {
                            "source": "#main/carveme/carveme_gem",
                            "id": "#main/memote_report_snapshot/GEM"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_report_snapshot/report_snapshot"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_report_snapshot/skip_test_find_metabolites_consumed_with_closed_bounds"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_report_snapshot/skip_test_find_metabolites_not_consumed_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_report_snapshot/skip_test_find_metabolites_not_produced_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_report_snapshot/skip_test_find_metabolites_produced_with_closed_bounds"
                        }
                    ],
                    "out": [
                        "#main/memote_report_snapshot/report_html"
                    ],
                    "id": "#main/memote_report_snapshot"
                },
                {
                    "label": "MEMOTE report snapshot",
                    "doc": "MEMOTE run analsis",
                    "run": "#memote.cwl",
                    "scatter": [
                        "#main/memote_run/GEM"
                    ],
                    "in": [
                        {
                            "source": "#main/carveme/carveme_gem",
                            "id": "#main/memote_run/GEM"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_run/run"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_run/skip_test_find_metabolites_consumed_with_closed_bounds"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_run/skip_test_find_metabolites_not_consumed_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_run/skip_test_find_metabolites_not_produced_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#main/memote_run/skip_test_find_metabolites_produced_with_closed_bounds"
                        }
                    ],
                    "out": [
                        "#main/memote_run/run_json"
                    ],
                    "id": "#main/memote_run"
                },
                {
                    "label": "prodigal",
                    "doc": "prodigal gene/protein prediction",
                    "run": "#prodigal.cwl",
                    "scatter": [
                        "#main/prodigal/input_fasta"
                    ],
                    "in": [
                        {
                            "source": "#main/bins",
                            "id": "#main/prodigal/input_fasta"
                        },
                        {
                            "default": true,
                            "id": "#main/prodigal/single_mode"
                        }
                    ],
                    "out": [
                        "#main/prodigal/predicted_proteins_faa"
                    ],
                    "id": "#main/prodigal"
                },
                {
                    "doc": "Preparation of workflow output files to a specific output folder",
                    "label": "Prodigal proteins to folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Prodigal_proteins",
                            "id": "#main/prodigal_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/compress_prodigal/outfile"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/prodigal_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/prodigal_files_to_folder/results"
                    ],
                    "id": "#main/prodigal_files_to_folder"
                },
                {
                    "label": "SMETANA",
                    "doc": "Species METabolic interaction ANAlysis",
                    "when": "$(inputs.run_smetana)",
                    "run": "#smetana.cwl",
                    "in": [
                        {
                            "source": "#main/carveme/carveme_gem",
                            "id": "#main/smetana/GEM"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/smetana/identifier"
                        },
                        {
                            "source": "#main/run_smetana",
                            "id": "#main/smetana/run_smetana"
                        },
                        {
                            "source": "#main/solver",
                            "id": "#main/smetana/solver"
                        }
                    ],
                    "out": [
                        "#main/smetana/detailed_output_tsv"
                    ],
                    "id": "#main/smetana"
                }
            ],
            "id": "#main",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2022-06-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        }
    ],
    "cwlVersion": "v1.2",
    "$namespaces": {
        "s": "https://schema.org/"
    }
}
