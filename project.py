from pyfiglet import Figlet
from bases import base_data
from residues import residue_data
from codons import codon_table


figlet = Figlet()


class DNA():
    valid_bases = "cgat"

    def __init__(self, dna_chain):
        for x in dna_chain:
            if x not in DNA.valid_bases:
                raise ValueError("Invalid DNA Base")
        self.chain = dna_chain
        self.length = len(dna_chain)
    # Initial State

    def content(self):
        count = 0
        total = self.length
        for x in self.chain:
            if x in "cg":
                count += 1
        content = count / total * 100
        return content
    # Content Method: Calculates CG Content
    
    def weight(self):
        weight = 0
        for x in self.chain:
            weight += base_data[x]["Molecular Weight"]
        return weight
    # Weight Method: Calcultes Molecular Weight

    def convert_1(self):
        print("\nRNA Chain: ", end="")
        print(self.chain.replace("t", "u"), end="")
    # Convert Method 1: Converts DNA Sequence To RNA Complement

    def convert_2(self):
        rna = RNA(self.chain.replace("t", "u"))
        print("")
        rna.convert()
    # Convert Method 2: Converts DNA Sequence To Peptide Chain
# DNA Class: For Analysing Or Converting DNA Sequences


class RNA():
    valid_bases = "cgau"

    def __init__(self, rna_chain):
        for x in rna_chain:
            if x not in RNA.valid_bases:
                raise ValueError("Invalid RNA Base")
        self.chain = rna_chain
        self.length = len(rna_chain)
    # Initial State

    def content(self):
        count = 0
        total = self.length
        for x in self.chain:
            if x in "cg":
                count += 1
        content = count / total * 100
        return content
    # Content Method: Calculates CG Content

    def weight(self):
        weight = 0
        for x in self.chain:
            weight += base_data[x]["Molecular Weight"]
        return weight
    # Weight Method: Calcultes Molecular Weight

    def convert(self):
        if self.chain.startswith("aug"):
            print("Peptide Chain: ", end="")
            for x in range(0, len(self.chain), 3):
                codon = self.chain[x:x+3]
                for k, v in codon_table.items():
                    if k == codon:
                        if v == "stop":
                            return
                        else:
                            print(f"{v}", end="")              
        else:
            raise ValueError("Invalid RNA Sequence")
    # Convert Method: Converts RNA Sequence To Peptide Chain
# RNA Class: For Analysing Or Converting RNA Sequences


class Peptide():
    valid_residues = "acdefghiklmnpqrstvwy"

    def __init__(self, peptide_chain):
        for x in peptide_chain:
            if x not in Peptide.valid_residues:
                raise ValueError("Invalid Residue")
        self.chain = peptide_chain
        self.length = len(peptide_chain)
    # Initial State

    def weight(self):
        weight = 0
        for x in self.chain:
            weight += residue_data[x]["Molecular Weight"]
        return weight
    # Weight Method: Calcultes Molecular Weight
# Peptide Class: For Analysing Peptides


class Residue():
    def __init__(self, residue_choice):
        count = 0
        for value in residue_data.values():
            if value["Name"] == residue_choice:
                count += 1
        if count == 1:
                self.choice = residue_choice
        else:
            raise ValueError("Invalid Base")
    # Initial State   
    
    def name(self):
        name = self.choice
        return name
    # Name Method: Returns Residue Name

    def code(self):
        for k, v in residue_data.items():
            if v["Name"] == self.choice:
                code = k.capitalize()
        return code
    # Code Method: Returns Residue 1-Letter Code

    def hydrophobicity(self):
        for value in residue_data.values():
            if value["Name"] == self.choice:
                hydrophobicity = value["Estimated Average Hydrophobicity"]
        return hydrophobicity
    # Hydrophobicity Method: Returns Residue Hydrophobicity

    def pi(self):
        for value in residue_data.values():
            if value["Name"] == self.choice:
                pi = value["pI"]
        return pi
    # PI Method: Returns Residue PI Value

    def weight(self):
        for value in residue_data.values():
            if value["Name"] == self.choice:
                weight = value["Molecular Weight"]
        return weight
    # Weight Method: Returns Molecular Weight
# Residue Class: For Analysing Individual Residues


class Base():
    def __init__(self, base_choice):
        count = 0
        for value in base_data.values():
            if value["Name"] == base_choice:
                count += 1
        if count == 1:
                self.choice = base_choice
        else:
            raise ValueError("Invalid Base")
    # Initial State

    def name(self):
        name = self.choice
        return name
    # Name Method: Returns Base Name

    def code(self):
        for value in base_data.values():
            if value["Name"] == self.choice:
                code = value["1-Letter Code"]
        return code
    # Code Method: Returns Base 1-Letter Code

    def type(self):
        for value in base_data.values():
            if value["Name"] == self.choice:
                type = value["Base Type"]
        return type
    # Type Method: Returns Base Type

    def complement(self):
        for value in base_data.values():
            if value["Name"] == self.choice:
                complement = value["Complement"]
        return complement
    # Complement Method: Returns Complement Base

    def hydrogen_bonds(self):
        for value in base_data.values():
            if value["Name"] == self.choice:
                hydrogen_bonds = value["Hydrogen Bonds"]
        return hydrogen_bonds
    # Hydrogen Bonds Method: Returns Number Of Hydrogen Bonds

    def weight(self):
        for value in base_data.values():
            if value["Name"] == self.choice:
                weight = value["Molecular Weight"]
        return weight
    # Weight Method: Returns Molecular Weight
# Base Class: For Analysing Individual Bases

def retry():
    retry = input("Press Enter To Continue: ")
    if retry != "Hello There!":
        main()
    else:
        print("\nGeneral Kenobi!")
    # Continue

def main():
    try:
        print(figlet.renderText("Peptide Analyser II"))
        choice = input(
            "'1': DNA Analyser And Converter"
            "\n'2': RNA Analyser And Converter"
            "\n'3': Peptide Analyser"
            "\n'4': Residue Analyser"
            "\n'5': Base Analyser" 
            "\n'6': DNA To RNA Converter" 
            "\n'7': DNA To Peptide Converter"
            "\n'8': RNA To Peptide Converter"
            "\n\nYour Choice: "
            ).strip().lower()
        # User Input Determines Pathway

        if choice == "1":
            dna_chain = input("\nEnter DNA Coding Sequence: ").strip().lower()
            dna = DNA(dna_chain)
            print(
            f"\nChain Length: {dna.length}"
            f"\nMolecular Weight: {dna.weight():.2f} Da"
            f"\nGC Content: {dna.content():.2f}%\n"
            )
        # Provides Information On DNA Chain Input

        elif choice == "2":
            rna_chain = input("\nEnter mRNA Sequence: ").strip().lower()
            rna = RNA(rna_chain)
            print(
            f"\nChain Length: {rna.length}"
            f"\nGC Content: {rna.content():.2f}%"
            f"\nMolecular Weight: {rna.weight():.2f} Da\n"
            )
        # Provides Information On RNA Chain Input   

        elif choice == "3":
            peptide_chain = input("\nEnter Peptide Chain: ").strip().lower()
            peptide = Peptide(peptide_chain)
            print(
            f"\nChain Length: {peptide.length}"
            f"\nMolecular Weight: {peptide.weight():.2f} Da\n"
            )
        # Provides Information On Peptide Chain Input

        elif choice == "4":
            residue_choice = input("\nEnter Residue: ").strip().lower()
            residue_choice = residue_choice.capitalize()
            residue = Residue(residue_choice)
            print(
            f"\nResidue: {residue.name()}"
            f"\n1-Letter Code: {residue.code()}"
            f"\nHydrophobicity: {residue.hydrophobicity()}"
            f"\npI: {residue.pi()}"
            f"\nMolecular Weight: {residue.weight():.2f} Da\n"
            )
        # Provides Information On Residue Input

        elif choice == "5":
            base_choice = input("\nEnter Base: ").strip().lower().capitalize()
            base_choice = base_choice.capitalize()
            base = Base(base_choice)
            print(
            f"\nBase: {base.name()}"
            f"\n1-Letter Code: {base.code()}"
            f"\nBase Type: {base.type()}"
            f"\nComplement Base: {base.complement()}"
            f"\nHydrogen Bonds: {base.hydrogen_bonds()}"
            f"\nMolecular Weight: {base.weight():.2f} Da\n"
            )
        # Provides Information On Base Input

        elif choice == "6":
            dna_chain = input("\nEnter DNA Coding Sequence: ").strip().lower()
            result = DNA(dna_chain)
            result.convert_1()
            print("\n")
        # Converts DNA Chain To RNA Chain
        
        elif choice == "7":
            dna_chain = input("\nEnter DNA Coding Sequence (Must Begin With A Valid Start Codon): ").strip().lower()
            result = DNA(dna_chain)
            result.convert_2()
            print("\n")
        # Converts DNA Chain To Peptide Chain

        elif choice == "8":
            rna_chain = input("\nEnter mRNA Sequence (Must Begin With A Valid Start Codon): ").strip().lower()
            print("")
            result = RNA(rna_chain)
            result.convert()
            print("\n")
        # Converts RNA Chain To Peptide Chain

        else:
            raise ValueError("Invalid Choice")
        # Raises Value Error If Non-Existant Choice Is Selected
        
    except:
        print("\nInput Error\n")
        # Catches Errors Neatly 

    input("\nClick Enter To Continue: ")
    main()
    # To Continue
# Main Function: Carries Out The Process


if __name__ == "__main__":     
    main()
# Runs Main If The Correct File Is Open

