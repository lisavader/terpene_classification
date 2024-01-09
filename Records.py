import re
import json
from ete3 import NCBITaxa
from types import SimpleNamespace

#query lists for enzyme classification
enzyme_types = {
    "terpenes": {
        "monoterpene" : ("monoterpene","iridoid","methylisoborneol","cineole"),
        "sesquiterpene" : ("sesquiterpene","sesquiterpenoid","trichodiene","trichodience","protoilludene","aristolochene","guaiadiene","germacrene"
                           ,"farnesene","viridiflorene","cuprenene","muurolene","longiborneol","presilphiperfolanol","eremophilene","cadinene"
                           ,"santalene","isozizaene","pentalenene","pristinol","isoafricanol","muurolol","cubebol","caryolanol","selinadiene"
                           ,"linalool","eudesmol","amorphene","corvol ether","germacradienol","caryophyllene","avermitilol"),
        "diterpene" : ("diterpene","variediene","fusicoccadiene","conidiogenone","araneosene","kaurene","phomopsene","dolastadiene","aphidicolanol"
                       ,"paxilline","BcPAX1","copalyl diphosphate","gibberellin","phyllocladanol","labdatriene","cyclooctatenol","tsukubadiene",
                       "spiroviolene","pimaradiene","terpentetriene"),
        "sesterterpene" : ("sesterterpene","ophiobolin","preasperterpenoid","preaspterpenacid","stellatatriene","astellifadiene","quiannulatene"
                           ,"sesterbrasiliatriene","sesterfisherol","preterpestacin","aspergildiene","fusoxypene"),
        "triterpene" : ("triterpene","macrophomene","talaropentaene"),
        "tetraterpene" : ("tetraterpene","lycopene")
    },
    "prenyltransferases": {
        "FPP" : ("FPP","farnesyl pyrophosphate","farnesyl diphosphate","dimethylallyltranstransferase","ERG20","geranyltranstransferase"),
        "GGPP" : ("GGPP","geranyl geranyl pyrophosphate","geranyl geranyl diphosphate","BTS1","farnesyltranstransferase"),
        "diapophytoene" : ("diapophytoene",),
        "phytoene" : ("phytoene",),
        "presqualene PP" : ("presqualene diphosphate",),
        "hydroxysqualene" : ("hydroxysqualene",),
        "squalene" : ("squalene","ERG9"),
        "HPP" : ("HPP","hexaprenyl pyrophosphate","hexaprenyl diphosphate","COQ1")
    }
}

class ProteinRecord:
    def __init__(self, database, header, accession, reviewed, protein_name, taxid):
        self.database = database
        self.header = header
        self.accession = accession
        self.reviewed = reviewed
        self.protein_name = protein_name
        self.taxid = taxid

        #initiate labels
        self.organism_name = ""
        self.organism_category = ""
        self.enzyme_type = set()      #can have multiple functions

    def assign_taxonomic_labels(self):
        ncbi = NCBITaxa()       #connector to ncbi database
        taxid_translator = ncbi.get_taxid_translator([self.taxid])
        self.organism_name = list(taxid_translator.values())[0]     #the first taxonomic name found for this taxid
        lineage = ncbi.get_lineage(self.taxid)
        if 5204 in lineage:
            self.organism_category = "Basidiomycota"
        elif 4890 in lineage:
            self.organism_category = "Ascomycota"
        elif 2 in lineage:
            self.organism_category = "Bacteria"
        else:
            self.organism_category = "Other"

    def assign_enzyme_type(self):
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
                    self.enzyme_type.add(terpene_type+" synthase")

        # Recognise prenyltransferases
        for pt_type, pt_names in enzyme_types["prenyltransferases"].items():
            for pt_name in pt_names:
                pattern = pt_name.replace(" ","[_-]?")
                match = re.search(pattern, self.protein_name, flags=re.IGNORECASE)
                if match:
                    self.enzyme_type.add(pt_type+" synthase")

        if not self.enzyme_type:
            self.enzyme_type.add("unknown")

        # Convert set to list (for compatibility with json)
        self.enzyme_type = list(self.enzyme_type)

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, indent=4)

    @staticmethod
    def from_json(json_file):
        with open(json_file, 'r') as handle:
            return json.load(handle, object_hook=lambda d: SimpleNamespace(**d))

class InterproRecord(ProteinRecord):
    def __init__(self, header):
        database = "interpro"
        header = header

        #split header into relevant fields
        fields = header.split("|")
        accession = fields[0]
        reviewed = fields[1] == "reviewed"
        protein_name = fields[2]
        taxid = int(fields[3].removeprefix("taxID:"))

        super().__init__(database, header, accession, reviewed, protein_name, taxid)
