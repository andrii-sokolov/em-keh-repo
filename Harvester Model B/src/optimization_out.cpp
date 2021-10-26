#include "optim.h"

std::string optimization::GNUPlot_string()
{
    return("plot \"ResonanceHB"+experimental_filename+"\" using 1:3 title \" Theory \" lt 1, \""+experimental_filename+"\" title \"Experiment\" lt 2");
}
