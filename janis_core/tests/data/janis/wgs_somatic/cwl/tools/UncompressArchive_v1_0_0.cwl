#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: UncompressArchive

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: ubuntu:latest

inputs:
- id: file
  label: file
  type: File
  inputBinding:
    position: 1
- id: stdout
  label: stdout
  doc: write on standard output, keep original files unchanged
  type: boolean
  default: true
  inputBinding:
    prefix: -c
- id: decompress
  label: decompress
  doc: decompress
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -d
- id: force
  label: force
  doc: force overwrite of output file and compress links
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -f
- id: keep
  label: keep
  doc: keep (don't delete) input files
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -k
- id: list
  label: list
  doc: list compressed file contents
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -l
- id: noName
  label: noName
  doc: do not save or restore the original name and time stamp
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -n
- id: name
  label: name
  doc: save or restore the original name and time stamp
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -N
- id: quiet
  label: quiet
  doc: suppress all warnings
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -q
- id: recursive
  label: recursive
  doc: operate recursively on directories
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -r
- id: suffix
  label: suffix
  doc: use suffix SUF on compressed files
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -s
- id: test
  label: test
  doc: test compressed file integrity
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -t
- id: fast
  label: fast
  doc: compress faster
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: '-1'
- id: best
  label: best
  doc: compress better
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: '-9'
- id: rsyncable
  label: rsyncable
  doc: Make rsync-friendly archive
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --rsyncable

outputs:
- id: out
  label: out
  type: stdout
stdout: _stdout
stderr: _stderr

baseCommand: gunzip
arguments: []
id: UncompressArchive
