# process the info part of vcf files from the 1000 genomes project
# stores the converted files in a directory named converted/ in current location

import os
import sys


def processFile(input_dir, input_fn, output_dir):
    print("Begin processing file: " + input_fn)
    
    output_fn = output_dir + input_fn[:input_fn.find('.')] + "_converted.vcf"
    output = open(output_fn, "w+")

    with open(os.path.join(input_dir, input_fn), "r") as f:
        for line in f:
            # parts: chromosome, position, ID, reference, alternate, quality, filter, info
            parts = line.split()

            # process info
            info = parts[7]
            info = info.split(';')


            info_dict = {}
            # split the info fields
            for field in info:
                field = field.rstrip().split('=')


                if len(field) == 2:
                    key, value = field

                    info_dict[key] = value


            # process the consequences portion of info
            if 'CSQ' in info_dict:
                csq = info_dict['CSQ'].split(',')

                for entry in csq:
                    entry = entry.split('|')

                    # only keep entry if it's the reported allele variant
                    if entry[0] == parts[4]:
                        # entry[0] = allele associated with entry's consequence
                        # parts[4] = alternate allele
                        
                        # use ONLY the first entry that matches the allele
                        info_dict["feature"] = entry[3]
                        info_dict["feature_type"] = entry[4]
                        break


            # we have all relevant fields
            # print output
            output_string = "{}\t{}\t{}\t{}\t{}\t{}".format(parts[0], parts[1], parts[2], parts[3], parts[4], info_dict["AF"])

            if "feature" in info_dict:
                output_string += "\t" + info_dict["feature"]
            else:
                output_string += "\t."


            if "feature_type" in info_dict:
                output_string += "\t" + info_dict["feature_type"]
            else:
                output_string += "\t."

            output.write(output_string + "\n")


    output.close()

    print("file completed")


def main():
    if len(sys.argv) < 2:
        print("One arguments are required:\n"
                "arg1: directory containing vcf files to convert."
                )
        exit()

    directory = sys.argv[1]

    output_dir = os.path.join(directory, "converted/")

    try:
        os.stat(output_dir)
    except:
        os.mkdir(output_dir)


    for fn in os.listdir(directory):
        if fn.endswith(".vcf"):
            processFile(directory, fn, output_dir)





if __name__ == "__main__":
    main()
