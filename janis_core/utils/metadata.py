import json
from _sha1 import sha1
from datetime import date


class Metadata(object):

    DATE_FORMAT = "%Y-%m-%d"

    def __init__(
        self,
        contributors=None,
        dateCreated=None,
        dateUpdated=None,
        institution=None,
        doi=None,
        citation=None,
        keywords=None,
        documentationUrl=None,
        documentation=None,
        short_documentation=None,
        version=None,
        sample_input_overrides: dict = None,
    ):
        """

        :param contributors:
        :param dateCreated:
        :param dateUpdated:
        :param institution:
        :param doi:
        :param citation:
        :param keywords:
        :param documentationUrl:
        :param documentation:
        :param short_documentation:
        :param version:
        """
        self.contributors = []
        if contributors:
            self.contributors = (
                contributors
                if isinstance(contributors, list)
                else contributors.split(",")
            )

        self.dateCreated = (
            dateCreated.strftime(self.DATE_FORMAT)
            if isinstance(dateCreated, date)
            else dateCreated
        )
        self.dateUpdated = (
            dateUpdated.strftime(self.DATE_FORMAT)
            if isinstance(dateUpdated, date)
            else dateUpdated
        )
        self.institution = institution
        self.doi = doi
        self.citation = citation
        self.keywords = keywords
        self.documentation = documentation
        self.short_documentation = short_documentation
        self.documentationUrl = documentationUrl
        self.version = version
        self.sample_input_overrides = sample_input_overrides

    def update(self, **kwargs):
        for k in kwargs:
            self.__setattr__(k, kwargs[k])
        return self

    def get_dict(self, object_to_checksum):

        checksum = sha1(
            json.dumps(object_to_checksum, sort_keys=True).encode("utf-8")
        ).hexdigest()
        # https://stackoverflow.com/q/5884066

        d = {k: v for k, v in vars(self).items() if v is not None}
        d["checksum"] = checksum
        d["dateGenerated"] = date.today().strftime(self.DATE_FORMAT)

        return d


class WorkflowMetadata(Metadata):
    def get_dict(self, object_to_checksum, ninputs=None, nsteps=None, noutputs=None):
        d = super(WorkflowMetadata, self).get_dict(object_to_checksum)

        if ninputs:
            d["numberOfInputs"] = ninputs
        if nsteps:
            d["numberOfSteps"] = nsteps
        if noutputs:
            d["numberOfOutputs"] = noutputs

        return d


class ToolMetadata(Metadata):
    pass
