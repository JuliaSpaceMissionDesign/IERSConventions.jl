# # [Loading IERS EOP data](@id tutorial_01_load)
# _This example was generated on DATEOFTODAY._

# When working with frames associated with the Earth, it is imperative to incorporate 
# Earth Orientation Parameters (EOP) data. The EOP data are required for the accurate 
# construction of  various Earth associated frames.

# ### Creating compatible EOP file

# To minimize dependencies on external sources, `IERSConventions` defines a 
# standardized format for EOP data. The expected format consists of a file with 
# the '.eop.dat' extension. This file should contain columns representing different 
# Earth orientation parameters. 

# In order to ensure that the provided EOP file adheres to this format, USNO standard files can 
# be processed using [`eop_generate_from_txt`](@ref) or [`eop_generate_from_csv`](@ref) 
# functions:

using IERSConventions

finals = download(
    "https://datacenter.iers.org/data/latestVersion/finals.all.iau2000.txt",
    "/tmp/finals2000A.txt"
)
filename = "/tmp/iau2000a"
eop_generate_from_txt(iers2010a, finals, filename)

# This will create a `iau2000a.eop.dat` file to be used later.

# ### Loading EOP data 

# Once a EOP file compatible with the reader is avaliable, the data could be 
# loaded in the environment for usage of all the methods within `IERSConventions` using 
# the [`eop_load_data!`](@ref) method, as follows:

eop_load_data!(iers2010a, filename * ".eop.dat")

# It is possible to override the data, e.g. re-use this functions as many times as needed.