#include "test_common.h"

static void test_status(void)
{
	printf("test_status\n");
	TEST_ASSERT(qaws_status_to_string(QAWS_STATUS_OK) != NULL, "status_to_string OK");
	TEST_ASSERT(qaws_status_to_string(QAWS_STATUS_INTERNAL_ERROR) != NULL, "status_to_string INTERNAL_ERROR");
}

int test_00_status_main(void)
{
	g_pass = 0; g_fail = 0;
	test_status();
	return g_fail > 0 ? 1 : 0;
}
