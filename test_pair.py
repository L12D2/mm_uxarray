from melodies_monet.driver import analysis
import traceback
import os

print("starting")

try:
    an = analysis()

    print("setting control")
    an.control = 'docs/examples/ungridded_support/control_camchem_se.yaml'

    print("reading control")
    an.read_control()

    print("control loaded")

    print("\nCONTROL DICT:")
    print(an.control_dict)

    print("\nstarting pair_data")
    an.pair_data()

    print("\npair_data finished")

    print("\nANALYSIS OBJECT ATTRIBUTES:")
    print(an.__dict__.keys())

    # Try common paired-data containers
    possible_attrs = [
        "paired",
        "pair",
        "models",
        "obs",
        "model",
        "df",
        "data"
    ]

    for attr in possible_attrs:
        if hasattr(an, attr):
            print(f"\nFOUND ATTRIBUTE: {attr}")
            obj = getattr(an, attr)

            print(type(obj))

            try:
                print(obj)
            except:
                print("Could not print object")

    outdir = "./output/airnow_camchemse"

    print(f"\nChecking output directory: {outdir}")

    if os.path.exists(outdir):
        print("Output directory exists")

        for root, dirs, files in os.walk(outdir):
            for f in files:
                print(os.path.join(root, f))
    else:
        print("Output directory does NOT exist")

    if hasattr(an, "plotting"):
        print("\nFound plotting method")
        print("Running plotting()")

        an.plotting()

        print("plotting finished")

    else:
        print("\nNo plotting() method found")

except Exception as e:
    print("\nERROR:")
    traceback.print_exc()
