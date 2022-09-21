nextflow.enable.dsl=2

// defaults  (set values in params.yaml)
params.greeting = null 
params.extra_greeting = null 

ch_greeting = Channel.of(params.greeting) 
ch_extra_greeting = Channel.of(params.extra_greeting) 
ch_false = Channel.of(true)

process CHECKFALSE {
    input:
    val input1
    cpus null

    output: 
    stdout emit: outstr

    script:
    def mytext = input1 ? 'its true!' : 'its false!'
    """ 
    echo $mytext
    """ 
}

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
    //results_ch = CONVERTTOUPPER(ch_greeting, ch_extra_greeting) 
    //results_ch.view{ it } 
    ch_out = CHECKFALSE(ch_false)
    ch_out.view{ it }
}


