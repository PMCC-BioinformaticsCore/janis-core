



from typing import Any, Optional

from janis_core.ingestion.galaxy.runtime.exceptions import AttributeNotSupportedError
from galaxy.tool_util.verify.interactor import ToolTestDescription
from .checks import ValidCheck
from .mapper import map_ttestout

from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput
)


"""
for each galaxy test, need to translate into a janis TTestCase. 
inputs are easy, but outputs are more complex. 
for each output, need to specify the TTestExpectedOutput.
this TTestExpectedOutput has a TTestPreprocessor, which obtains the value of the
"""

class TestFactory:
    def produce(self, gxtest: ToolTestDescription) -> TTestCase:
        return TTestCase(
            str(gxtest.name),
            self.prepare_inputs(gxtest),
            self.prepare_outputs(gxtest)
        )

    def prepare_inputs(self, gxtest: ToolTestDescription) -> dict[str, Any]:
        """maps the galaxy test def inputs to janis format. not many changes."""
        inputs: dict[str, Any] = {}
        for key, val in gxtest.inputs.items():
            if isinstance(val, list) and len(val) == 1:
                val = val[0]
            inputs[str(key)] = val
        return inputs
        
    def prepare_outputs(self, gxtest: ToolTestDescription) -> list[TTestExpectedOutput]:
        outputs: list[TTestExpectedOutput] = []
        for gxout in gxtest.outputs:
            checks = self.gather_checks(gxout)
            outputs += [map_ttestout(c) for c in checks]
        return outputs

    def gather_checks(self, gxout: dict[str, Any]) -> list[ValidCheck]:
        checks: list[ValidCheck] = []
        checks += self.gather_diff_checks(gxout)
        checks += self.gather_re_match_checks(gxout)
        checks += self.gather_sim_size_checks(gxout)
        checks += self.gather_file_fingerprint_checks(gxout)
        checks += self.gather_file_assert_checks(gxout)
        checks += self.gather_file_contains_checks(gxout)
        return checks

    def gather_diff_checks(self, gxout: dict[str, Any]) -> list[ValidCheck]:
        if gxout['attributes']['compare'] == 'diff' and gxout['value']:
            return [self.create_check(gxout['name'], 'lines_diff', gxout['attributes']['lines_diff'], reffile=gxout['value'])]
        return []
    
    def gather_re_match_checks(self, gxout: dict[str, Any]) -> list[ValidCheck]:
        if gxout['attributes']['compare'] in ['re_match', 're_match_multiline']:
            raise AttributeNotSupportedError('re_match and re_match_multiline not yet supported')
        return []
    
    def gather_sim_size_checks(self, gxout: dict[str, Any]) -> list[ValidCheck]:
        if gxout['attributes']['compare'] == 'sim_size':
            fields = ['delta', 'delta_frac']
            return self.create_attribute_checks(gxout, valid_fields=fields)
        return []
        
    def gather_file_fingerprint_checks(self, gxout: dict[str, Any]) -> list[ValidCheck]:
        # file fingerprints
        # TODO add ftype checks!
        #fields = ['ftype', 'md5', 'checksum']
        fields = ['md5', 'checksum']

        return self.create_attribute_checks(gxout, valid_fields=fields)

    def gather_file_assert_checks(self, gxout: dict[str, Any]) -> list[ValidCheck]:
        checks: list[ValidCheck] = []
        if gxout['attributes']['assert_list']:
            for assertcheck in gxout['attributes']['assert_list']:
                checks.append(self.create_check(gxout['name'], assertcheck['tag'], assertcheck['attributes']))
        return checks
        
    def gather_file_contains_checks(self, gxout: dict[str, Any]) -> list[ValidCheck]:
        checks: list[ValidCheck] = []
        if gxout['attributes']['compare'] == 'contains' and gxout['value']:
            checks.append(self.create_check(gxout['name'], 'has_text', None, reffile=gxout['value']))
        return checks

    def create_attribute_checks(self, gxout: dict[str, Any], valid_fields: Optional[list[str]] = None) -> list[ValidCheck]:
        """creates an attribute checks based on the values of valid_valid_fields"""
        if not valid_fields:
            valid_fields = []
        return [ self.create_check(gxout['name'], attname, attval) 
                 for attname, attval in gxout['attributes'].items() 
                 if attname in valid_fields and attval]

    def create_check(self, name: str, ctype: str, value: Optional[Any], reffile: Optional[str]=None) -> ValidCheck:
        return ValidCheck(name, ctype, value, reffile)

    def create_check2(self, name: str, ctype: str, value: Optional[Any], reffile: Optional[str]=None) -> ValidCheck:
        if value:
            if not isinstance(value, dict):
                test_values: dict[str, Any] = {'value': value}
            else:
                test_values = value
        else:
            test_values: dict[str, Any] = {}
        return ValidCheck(name, ctype, test_values, reffile)


        





    #     # get the right map object 
    #     ttestouts = []
    #     for out in gxtest.outputs:
    #         ttestouts += self.init_file_diff_outputs(out)
    #         ttestouts += self.init_file_check_outputs(out)
    #     return ttestouts

    # def init_file_diff_outputs(self, gx_test_out: dict[str, Any]) -> list[TTestExpectedOutput]:
    #     outs = []
    #     if gx_test_out['value']:
    #         self.init_diff_outputs(gx_test_out)
    #     return outs

    # def init_file_check_outputs(self, gx_test_out: dict[str, Any]) -> list[TTestExpectedOutput]:
    #     outs = []
    #     outs += self.init_fingerprint_outputs(gx_test_out)
    #     outs += self.init_rematch_outputs(gx_test_out)
    #     outs += self.init_simsize_outputs(gx_test_out)
    #     outs += self.init_assert_outputs(gx_test_out)
    #     return outs

    # def init_fingerprint_outputs(self, gx_test_out: dict[str, Any]) -> list[TTestExpectedOutput]:
    #     pass

    # def init_fingerprint_outputs(self, gx_test_out: dict[str, Any]) -> list[TTestExpectedOutput]:
    #     if gx_test_out['attributes']['compare'] in ['re_match', 're_match_multiline']:
    #         raise AttributeNotSupportedError('re_match and re_match_multiline output checks not supported')

    # def init_simsize_outputs(self, gx_test_out: dict[str, Any]) -> list[TTestExpectedOutput]:
    #     outchecks = ['']
    #     outs = []
    #     for check in outchecks:
    #         outs += strategy_map[check](gx_test_out)
    #     return outs

    # def init_assert_outputs(self, gx_test_out: dict[str, Any]) -> list[TTestExpectedOutput]:
    #      if gx_test_out['attributes']['assert_list']:
    #         for out_check in gx_test_out['attributes']['assert_list']:
    #             self.map_assert_output()





"""
COMPARE FILE NEEDED -> may need decompression if decompress="true"
compare="diff":
- lines_diff
- ftype

compare="sim_size":
do filesize check with the following info
- delta
maximum allowed absolute size difference (in bytes) between the data set that is generated in the test and the file in test-data/ that is referenced by the file attribute. Default value is 10000 bytes. Can be combined with delta_frac.
- delta_frac
If compare is set to sim_size, this is the maximum allowed relative size difference between the data set that is generated in the test and the file in test-data/ that is referenced by the file attribute. A value of 0.1 means that the file that is generated in the test can differ by at most 10% of the file in test-data. The default is not to check for relative size difference. Can be combined with delta.


FILE NOT NEEDED
md5
checksum
assert_list has elems



"""





    # def init_output(self):
    #     pass
    #     for output in gxtest.outputs:
    #         self.
    #     if gxtest.outputs
    #     #preprocessor_map = self.
    #     compare_method = gxtest.output_data['compare']
    #     return TTestExpectedOutput(
    #         gxtest.name,
    #         preprocessor_map[compare_method],
    #     )
        
        



