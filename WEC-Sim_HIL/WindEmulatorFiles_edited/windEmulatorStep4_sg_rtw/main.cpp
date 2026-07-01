/* Main generated for Simulink Real-Time model windEmulatorStep4 */
#include <ModelInfo.hpp>
#include <utilities.hpp>
#include "rte_windEmulatorStep4_parameters.h"
#include "windEmulatorStep4.h"

/* Task wrapper function definitions */
void windEmulatorStep4_Task1(void)
{ 
    windEmulatorStep4_step();
} 
/* Task descriptors */
slrealtime::TaskInfo task_1( 0u, std::bind(windEmulatorStep4_Task1), slrealtime::TaskInfo::PERIODIC, 0.004, 0, 40);

/* Executable base address for XCP */
#ifdef __linux__
extern char __executable_start;
static uintptr_t const base_address = reinterpret_cast<uintptr_t>(&__executable_start);
#else
/* Set 0 as placeholder, to be parsed later from /proc filesystem */
static uintptr_t const base_address = 0;
#endif

/* Model descriptor */
slrealtime::ModelInfo windEmulatorStep4_Info =
{
    "windEmulatorStep4",
    windEmulatorStep4_initialize,
    windEmulatorStep4_terminate,
    []()->char const*& { return windEmulatorStep4_M->errorStatus; },
    []()->unsigned char& { return windEmulatorStep4_M->Timing.stopRequestedFlag; },
    { task_1 },
    slrealtime::getSegmentVector()
};

int main(int argc, char *argv[]) {
    slrealtime::BaseAddress::set(base_address);
    return slrealtime::runModel(argc, argv, windEmulatorStep4_Info);
}
