// Register ApprovalTests with Catch2 v3.
// This must be compiled as a separate translation unit.
#define APPROVALS_CATCH2_V3
#include "ApprovalTests.hpp"
#include "ApprovalTests/integrations/catch/Catch2v3Approvals.h"

auto directoryDisposer = ApprovalTests::Approvals::useApprovalsSubdirectory("approval_tests");
