#include "tbkp_instance.h"
#include <stdio.h>

int main() {
    TBKPInstance* instance = tbkp_instance_read("../data/example.txt");
    tbkp_instance_print(instance);
    tbkp_instance_free(&instance);

    return 0;
}
