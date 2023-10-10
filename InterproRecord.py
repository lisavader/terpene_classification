import re
import json
from ete3 import NCBITaxa
from types import SimpleNamespace

#query lists for enzyme classification
enzyme_types = {
    "terpenes": {
        "monoterpene" : ("monoterpene","iridoid"),
        "sesquiterpene" : ("sesquiterpene","sesquiterpenoid","trichodiene","trichodience","protoilludene","aristolochene","guaiadiene","germacrene"
                           ,"farnesene","viridiflorene","cuprenene","muurolene","longiborneol","presilphiperfolanol","eremophilene","cadinene"
                           ,"santalene","isozizaene","pentalenene","pristinol"),
        "diterpene" : ("diterpene","variediene","fusicoccadiene","conidiogenone","araneosene","kaurene","phomopsene","dolastadiene","aphidicolanol"
                       ,"paxilline","BcPAX1","copalyl diphosphate","gibberellin","phyllocladanol"),
        "sesterterpene" : ("sesterterpene","ophiobolin","preasperterpenoid","preaspterpenacid","stellatatriene","astellifadiene","quiannulatene"
                           ,"sesterbrasiliatriene","sesterfisherol","preterpestacin","aspergildiene","fusoxypene"),
        "triterpene" : ("triterpene","macrophomene","talaropentaene"),
        "tetraterpene" : ("tetraterpene","lycopene")
    },
    "prenyltransferases": {
        "FPP" : ("FPP","farnesyl pyrophosphate","farnesyl diphosphate","dimethylallyltranstransferase","ERG20","geranyltranstransferase"),
        "GGPP" : ("GGPP","geranyl geranyl pyrophosphate","geranyl geranyl diphosphate","BTS1","farnesyltranstransferase"),
        "phytoene" : ("phytoene",),
        "squalene" : ("squalene","ERG9"),
        "HPP" : ("HPP","hexaprenyl pyrophosphate","hexaprenyl diphosphate","COQ1")
    }
}

class InterproRecord:
    def __init__(self, header, seq):
        self.header = header
        self.seq = seq

        #split header into relevant fields
        fields = self.header.split("|")
        self.accession = fields[0]
        self.review_status = fields[1]
        self.protein_name = fields[2]
        self.taxid = int(fields[3].removeprefix("taxID:"))

        #initiate labels
        self.organism_name = ""
        self.organism_division = ""
        self.enzyme_class = set()      #can have multiple classes and subclasses
        self.enzyme_subclass = set()

    def assign_taxonomic_labels(self):
        ncbi = NCBITaxa()       #connector to ncbi database
        taxid_translator = ncbi.get_taxid_translator([self.taxid])
        self.organism_name = list(taxid_translator.values())[0]     #the first taxonomic name found for this taxid
        lineage = ncbi.get_lineage(self.taxid)
        if 5204 in lineage:
            self.organism_division = "Basidiomycota"
        elif 4890 in lineage:
            self.organism_division = "Ascomycota"
        else:
            self.organism_division = "Other"

    def assign_enzyme_labels(self):
        # Recognise terpene synthases
        for terpene_type, terpene_names in enzyme_types["terpenes"].items():
            for terpene_name in terpene_names:
                #spaces between words are replaced by an underscore or hyphen
                pattern = terpene_name.replace(" ","[_-]?")
                #compound suffixes -diene, -ene and -ol can be preceded by structural denotions between hyphens, e.g. dolasta-1(15),8-diene
                for suffix in ["diene","ene","ol"]:
                    if pattern.endswith(suffix):
                        pattern = re.sub(suffix, "(-.*-)?"+suffix, pattern)
                match = re.search(pattern, self.protein_name, flags=re.IGNORECASE)
                if match:
                    self.enzyme_class.add("terpene synthase")
                    self.enzyme_subclass.add(terpene_type+" synthase")
                else:       # if no specific terpene name is found, search for generic descriptors of a terpene synthase
                    for pattern in ["terpene", "terpenoid"]:
                        match = re.search(pattern, self.protein_name, flags=re.IGNORECASE)
                        if match:
                            self.enzyme_class.add("terpene synthase")

        # Recognise prenyltransferases
        for pt_type, pt_names in enzyme_types["prenyltransferases"].items():
            for pt_name in pt_names:
                pattern = pt_name.replace(" ","[_-]?")
                match = re.search(pattern, self.protein_name, flags=re.IGNORECASE)
                if match:
                    self.enzyme_class.add("prenyltransferase")
                    self.enzyme_subclass.add(pt_type+" synthase")
                else:
                    for pattern in ["prenyl transferase","polyprenyl","polyprenyl diphosphate"]:
                        match = re.search(pattern, self.protein_name, flags=re.IGNORECASE)
                        if match:
                            self.enzyme_class.add("prenyltransferase")

        if not self.enzyme_class:
            self.enzyme_class.add("unknown")
        if not self.enzyme_subclass:
            self.enzyme_subclass.add("unknown")

        # Convert set to list (for compatibility with json)
        self.enzyme_class = list(self.enzyme_class)
        self.enzyme_subclass = list(self.enzyme_subclass)

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, indent=4)

    @staticmethod
    def from_json(json_file):
        with open(json_file, 'r') as handle:
            return json.load(handle, object_hook=lambda d: SimpleNamespace(**d))
