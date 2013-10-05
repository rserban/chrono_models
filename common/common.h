#ifndef HAMMADMODELS_COMMON_H
#define HAMMADMODELS_COMMON_H
#include <vector>
#include "lcp/ChlcpVariablesGeneric.h"
#include "lcp/ChlcpVariablesBody.h"
#include "lcp/ChlcpConstraintTwoGeneric.h"
#include "lcp/ChlcpSystemDescriptor.h"
#include "lcp/ChlcpIterativeSOR.h"
#include "lcp/ChlcpIterativePMINRES.h"
#include "lcp/ChlcpIterativeAPGD.h"
#include "lcp/ChlcpIterativeBB.h"
#include "lcp/ChlcpSimplexSolver.h"
#include "core/ChlinearAlgebra.h"
#include "physics/ChsystemOpenMP.h"
#include "assets/ChsphereShape.h"
#include "physics/Chapidll.h"
#include "physics/Chsystem.h"
#include "lcp/ChlcpIterativeMINRES.h"
#include "core/ChrealtimeStep.h"
#include "lcp/ChlcpIterativeSOR.h"
#include "lcp/ChlcpIterativeJacobi.h"
#include "collision/ChcModelBullet.h"
#include "collision/ChcCollisionSystemBullet.h"
#include "physics/ChcontactContainer.h"
#include "ChopenGL.h"
#include "ChsystemParallel.h"
#include "ChlcpSystemDescriptorParallel.h"
#include "unit_POSTPROCESS/ChmitsubaRender.h"
//#include "unit_GPU/CHbodyFluid.h"
using namespace chrono;
using namespace postprocess;
using namespace geometry;
using namespace std;


#define PI 3.14159265359


#endif
