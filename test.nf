nextflow.enable.dsl=2

// defaults  (set values in params.yaml)
params.greeting = null 
params.extra_greeting = null 

ch_greeting = Channel.of(params.greeting) 
ch_extra_greeting = Channel.of(params.extra_greeting) 


process CONVERTTOUPPER { 
    input:  
    val text1
    val text2

    output: 
    stdout 

    script:
    def final_text = text2 ? text1 + ' ' + text2 : text1
    """ 
    echo "$final_text" | tr '[a-z]' '[A-Z]'  
    """ 
} 

workflow { 
    results_ch = CONVERTTOUPPER(ch_greeting, ch_extra_greeting) 
    results_ch.view{ it } 
}


