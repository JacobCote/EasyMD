
from openbabel import openbabel

def get_charges(input_file, output_file, file_type_input, file_type_output):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(file_type_input, file_type_output)

    #
    # 8JF
    #
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_file)

    # Calcul des charges MMFF94 (compatible avec ledock)
    charge_model = openbabel.OBChargeModel.FindType("MMFF94")
    charge_model.ComputeCharges(mol)

    obConversion.WriteFile(mol, output_file)

    print(f"Charges for ligand calculated and saved in {output_file}")


import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        prog='openbabel charge',
                        description='Calculates charges using openbabel and MMFF94 force field',
                        epilog='Enjoy the program! :)',)

    parser.add_argument('-i', '--input', type=str, required=True, help='Input file')
    parser.add_argument('-f', '--file_type', type=str, required=True, help='Input file format')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')


    args = parser.parse_args()
    get_charges()

    print(f"Calculating charges for {args.input} with MMFF94 and saving in {args.output}")



