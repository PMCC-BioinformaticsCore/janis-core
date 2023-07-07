
cwlVersion: v1.0
class: CommandLineTool

requirements:
    InlineJavascriptRequirement: {}
    ResourceRequirement:
      coresMin: 8
      ramMin: 3800

baseCommand: [ echo ]

inputs:
    # mandatory
    inSecondary1:
        type: File
        inputBinding:
            position: 6
            prefix: "-r"
        secondaryFiles: [".bai"]
    inSecondary2:
        type: File
        secondaryFiles: |
            ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                    return {"location": dpath+self.nameroot+".bai", "class": "File"}
                } 
                else{
                    return {"location": dpath+self.basename+".crai", "class": "File"}
                }
            } 
        doc: "tumor BAM or CRAM"
    inSecondaryArr1:
        type: File[]
        inputBinding:
            position: 6
            prefix: "-r"
        secondaryFiles: [".bai"]
    inSecondaryArr2:
        type: File[]
        secondaryFiles: |
            ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                    return {"location": dpath+self.nameroot+".bai", "class": "File"}
                } 
                else{
                    return {"location": dpath+self.basename+".crai", "class": "File"}
                }
            } 
        doc: "tumor BAM or CRAM"
    
    # optional
    inSecondaryOpt1:
        type: File?
        inputBinding:
            position: 6
            prefix: "-r"
        secondaryFiles: [".bai"]
    inSecondaryOpt2:
        type: File?
        secondaryFiles: |
            ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                    return {"location": dpath+self.nameroot+".bai", "class": "File"}
                } 
                else{
                    return {"location": dpath+self.basename+".crai", "class": "File"}
                }
            }
        doc: "tumor BAM or CRAM"
    inSecondaryArrOpt1:
        type: File[]?
        inputBinding:
            position: 6
            prefix: "-r"
        secondaryFiles: [".bai"]
    inSecondaryArrOpt2:
        type: File[]?
        secondaryFiles: |
            ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                    return {"location": dpath+self.nameroot+".bai", "class": "File"}
                } 
                else{
                    return {"location": dpath+self.basename+".crai", "class": "File"}
                }
            }
        doc: "tumor BAM or CRAM"

outputs: []