
#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
// #include "../src/imtiaz_2002d_noTstart_COR.hpp"
#include "../src/CellICCBioPhy.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
// #include "../src/DummyDerivedCa.hpp"
#include "Debug.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "BidomainProblem.hpp"
#include "ChasteEllipsoid.hpp"
#include "AbstractElement.hpp"

#include "FileComparison.hpp"
#include "FileFinder.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"

// for mechanics
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "QuadraticMesh.hpp"
#include "NonlinearElasticityTools.hpp"

#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"

#include "ExponentialMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
// cell factories
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"

// #include <GaussianQuadratureRule.hpp>
// #include <QuadraturePointsGroup.hpp>



using namespace std;

class ICCCellFactory : public AbstractCardiacCellFactory<3>
{

public:
  ICCCellFactory():AbstractCardiacCellFactory<3>()
  {
  }

  AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
  {
    // Cellimtiaz_2002d_noTstart_CORFromCellML* cell = new Cellimtiaz_2002d_noTstart_CORFromCellML(mpSolver, mpZeroStimulus);
    // cell->SetParameter("eta", 0.04);

    CellICCBioPhy* cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
    cell->SetParameter("V_excitation", -55); // excitation value (threshold) default: -55mV
    cell->SetParameter("live_time", 10000); // time of resting potential
    cell->SetParameter("ode_time_step", 0.1); // Set the same as defined in HeartHandler
    cell->SetParameter("IP3Par", 0.00069); //
    cell->SetParameter("t_start", 600000); // Set larger than total simulation time

    double x = pNode->rGetLocation()[0];
    double y = pNode->rGetLocation()[1];
    double z = pNode->rGetLocation()[2];
    // double r = 0.1; // set size of the radius
    // double scale = 1; // Set the same like for mesh.scale
    // ChastePoint<2> centre(0.1,0.45);
    // ChastePoint<2> radii (0.01,0.01);
    // ChasteEllipsoid<2> ellipseRegion(centre, radii);
    // ChastePoint<2> myPoint(x, y);

    ChastePoint<3> centre(0.39, -1.38, -3.13);
    ChastePoint<3> radii (0.05,0.05,0.05);
    ChasteEllipsoid<3> ellipseRegion(centre, radii);
    ChastePoint<3> myPoint(x, y, z);
    // double ca = cell->GetIntracellularCalciumConcentration();
    // std::cout << "Ca: " << ca << std::endl;
    if(ellipseRegion.DoesContain(myPoint))
    {
      cell->SetParameter("t_start", 0);
      std::cout<<"I'm inside region!" << std::endl;
    }
    return cell;
  }
};

class Test3DMonodomain : public CxxTest::TestSuite
{
public:
  void TestSimulation() //throw(Exception)
  {
    // Read electric mesh file
    TetrahedralMesh<3,3> mesh;
    std::string myFile = "scaffold_3elems_16_16_1.3";
    std::string meshFile = "/home/chaste/projects/mesh/3elems_16_16_1_v2/" + myFile;
    TrianglesMeshReader<3,3> mesh_reader(meshFile.c_str());
    mesh.ConstructFromMeshReader(mesh_reader);
    // mesh.Scale(20, 20, 20);

    // Read mechanics mesh f1ile
    QuadraticMesh<3> mechanics_mesh;
    std::string myFile_m = "scaffold_3elems_16_16_1.2";
    std::string meshFile_m = "/home/chaste/projects/mesh/3elems_16_16_1_v2/" + myFile_m;
    TrianglesMeshReader<3,3> mesh_reader_m(meshFile_m.c_str(),2);
    mechanics_mesh.ConstructFromMeshReader(mesh_reader_m);
    // mechanics_mesh.Scale(20, 20, 20);

    // TetrahedralMesh<3,3> mesh;
    // mesh.ConstructRegularSlabMesh(0.1/*stepsize*/, 0.2/*length*/, 0.5/*width*/, 0.1/*depth*/);
    //
    // QuadraticMesh<3> mechanics_mesh;
    // mechanics_mesh.ConstructRegularSlabMesh(0.1, 0.2, 0.5, 0.1 /*as above with a different stepsize*/);
    //
    //
    ///// fixed nodes are all nodes on top and bottom
    std::vector<unsigned> fixed_nodes;
    fixed_nodes  = {17,1636, 1649, 1598};
    // fixed_nodes = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 1, 0.0);

    // std::vector<unsigned> fixed_nodes_bottom  = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 1, 0.0);
    // std::vector<unsigned> fixed_nodes_top     = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 1, 0.5);
    //
    // for(unsigned int j = 0; j < fixed_nodes_top.size(); j++)
    // {
    //   fixed_nodes.push_back(fixed_nodes_top[j]);
    //   fixed_nodes.push_back(fixed_nodes_bottom[j]);
    // }

    // dev test
    std::cout << "fixed node size: " << fixed_nodes.size() << std::endl;
    for(unsigned int k = 0; k < fixed_nodes.size(); k++)
    std::cout << "fixed_nodes:" << fixed_nodes[k] << ' ' << std::endl;


    int counter = 0;
    // dev test
    for (unsigned i=0; i<mechanics_mesh.GetNumNodes()-1; i++)
    {
      bool bn = mechanics_mesh.GetNode(i)->IsBoundaryNode();
      double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
      double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
      double Z = mechanics_mesh.GetNode(i)->rGetLocation()[2];
      // std::cout <<"node id: " << i << ", X: " << X << " , Y: " << Y << " , Z: " << Z << std::endl;
      // int N = mechanics_mesh.GetNode(i)->GetIndex();
      // std::cout <<"node index: " << N << std::endl;
      if (bn){
        counter++;
         std::cout <<"I'm boundary node , node id: " << i << ", X: " << X << " , Y: " << Y << " , Z: " << Z << std::endl;
       }
    }
    std::cout <<"Number of boundary nodes: " << counter << std::endl;
    // Simulation settings
    HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
    HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
    HeartConfig::Instance()->SetCapacitance(2.5);
    HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12, 0.12));
    HeartConfig::Instance()->SetSimulationDuration(1000);  //ms.
    HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,100); // doesn't work here!

    HeartConfig::Instance()->SetPrintingTimeStep(100.0);
    double print_time_step = HeartConfig::Instance()->GetPrintingTimeStep();
    cout << "Print time step: " << print_time_step << endl;

    // Output visualization options
    HeartConfig::Instance()->SetVisualizeWithCmgui(false);
    HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
    HeartConfig::Instance()->SetVisualizeWithVtk(false);

    // Cell factory
    ICCCellFactory cell_factory;

    // Material law
    MooneyRivlinMaterialLaw<3> law(3.0, 3.0);
    ExponentialMaterialLaw<3> law2(1.5, 3.0); // First parameter is 'a', second 'b', in W=a*exp(b(I1-3))

    // Fibre directions
    // coming soon!

    // Cardiac ElectroMechanics problem definition
    ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
    problem_defn.SetContractionModel(NASH2004,0.1/*contraction model ODE timestep*/);
    // problem_defn.SetSolverType(IMPLICIT);
    // problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, false /*deformation affects cell models*/);
    // problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
    problem_defn.SetMaterialLaw(INCOMPRESSIBLE, &law2);
    problem_defn.SetZeroDisplacementNodes(fixed_nodes);
    problem_defn.SetMechanicsSolveTimestep(1.0);
    // problem_defn.SetVariableFibreSheetDirectionsFile(finder, true);
    problem_defn.SetSolveUsingSnes();

    CardiacElectroMechanicsProblem<3,1> problem(INCOMPRESSIBLE,
      MONODOMAIN,
      &mesh,
      &mechanics_mesh,
      &cell_factory,
      &problem_defn,
      "icc3d");

      c_vector<double,3> node_to_watch;
      node_to_watch(0) = 0.128477;
      node_to_watch(1) = -1.4124;
      node_to_watch(2) = -3.00469;
      problem.SetWatchedPosition(node_to_watch);
      problem.SetOutputDeformationGradientsAndStress(100.0);

      // solve
      // problem.Initialise();
      problem.Solve();


    }
  };
