# Converting CWL into WDL

This script parses a CWL command line tool, and converts it into a Janis CommandTool. It uses CWL-Utils to do the parsing,
and then a pretty simple map back into Janis. This is a proof of concept, with the intent to turn this script into something
more stable and useful (rather than just converting CWL -> WDL).

## Usage

This will try to guess the CWLVersion from the `cwlVersion` field on the command line tool. It only recognises

```
./fromcwl.py bwamem.cwl
```

This will print the WDL representation of the tool. 


### Requirements

```
pip install 'cwl-utils >= 0.6' 'ruamel.yaml >= 0.12.4, <= 0.16.5'
```


## Limitations

Started simple, and can add more as we go:

- Doesn't convert JavaScript expressions (even ones that Janis generates)
- Doesn't convert requirements (except Docker.pull)


 
