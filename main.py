import os
import shutil
from sys import stdout
from typing import Tuple

import openmm
from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box
from openff.toolkit import Molecule as openffMolecule
from openff.toolkit import Topology as openffTopology
from openff.units import unit as openff_unit
from openmm import app
from openmmforcefields.generators import SMIRNOFFTemplateGenerator


def prepare_openmm_system(
    host_guest_code: Tuple, ff: openmm.app.ForceField
) -> (openmm.System, openmm.app.Modeller):
    host_code, guest_code = host_guest_code
    host_mol = openffMolecule.from_file(
        f"{hosts_dir}/{host_code}/{host_code}.sdf", allow_undefined_stereo=True
    )
    guest_mol = openffMolecule.from_file(
        f"{hosts_dir}/{host_code}/{guest_code}/{guest_code}.sdf",
        allow_undefined_stereo=True,
    )

    smirnoff = SMIRNOFFTemplateGenerator(molecules=[host_mol, guest_mol])
    ff.registerTemplateGenerator(smirnoff.generator)

    for name in os.listdir(f"{hosts_dir}/{host_code}/{guest_code}"):
        if name[-5:] == "p.pdb":
            complex_topology = openffTopology.from_pdb(
                f"{hosts_dir}/{host_code}/{guest_code}/{name}",
                unique_molecules=[host_mol, guest_mol],
            )
            os.makedirs(f"working_data/{host_code}-{guest_code}", exist_ok=True)
            shutil.copy(
                f"{hosts_dir}/{host_code}/{guest_code}/{name}",
                f"working_data/{host_code}-{guest_code}/{name}",
            )

    for atom in list(complex_topology.molecules)[0].atoms:
        atom.metadata["residue_number"] = 1

    model = app.Modeller(
        complex_topology.to_openmm(ensure_unique_atom_names=False),
        complex_topology.get_positions().m_as(openff_unit.nanometer),
    )
    model.addSolvent(
        ff,
        padding=2.0 * openmm.unit.nanometers,
        ionicStrength=0.1 * openmm.unit.molar,
        positiveIon="Na+",
        negativeIon="Cl-",
    )
    system = ff.createSystem(model.topology, nonbondedMethod=app.PME)
    return system, model


def simulate(
    host_guest_code: Tuple, system: openmm.System, model: openmm.app.Modeller
) -> None:
    host_code, guest_code = host_guest_code
    working_dir = f"working_data/{host_code}-{guest_code}"

    timestep = 2 * openmm.unit.femtoseconds
    time_to_simulate = 10 * openmm.unit.nanosecond
    total_steps = int(time_to_simulate / timestep)

    integrator = openmm.LangevinMiddleIntegrator(
        298.15 * openmm.unit.kelvin, 1 / openmm.unit.picosecond, timestep
    )
    integrator.setRandomNumberSeed(12345678)
    platform = openmm.Platform.getPlatformByName("CUDA")
    simulation = app.Simulation(model.topology, system, integrator, platform)
    simulation.context.setPositions(model.positions)
    simulation.minimizeEnergy(
        0.1 * openmm.unit.kilojoules_per_mole / openmm.unit.nanometer
    )
    simulation.reporters.append(
        app.StateDataReporter(
            stdout,
            10000,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            totalSteps=total_steps,
            speed=True,
            progress=True,
        )
    )
    simulation.reporters.append(
        app.StateDataReporter(
            open(f"{working_dir}/production.csv", "w"),
            10000,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            speed=True,
        )
    )
    simulation.reporters.append(app.DCDReporter(f"{working_dir}/production.dcd", 10000))
    simulation.step(total_steps)


if __name__ == "__main__":
    host_guest_codes = []

    hosts_dir = "host-guest-benchmarks/taproom/systems"
    host_codes = [
        name
        for name in os.listdir(hosts_dir)
        if os.path.isdir(os.path.join(hosts_dir, name))
    ]

    for host_code in host_codes:
        guests_dir = f"{hosts_dir}/{host_code}"
        guest_codes = [
            name
            for name in os.listdir(guests_dir)
            if os.path.isdir(os.path.join(guests_dir, name))
        ]
        for guest_code in guest_codes:
            host_guest_code = (host_code, guest_code)
            host_guest_codes.append(host_guest_code)

    ff = app.ForceField("amber/tip3p_standard.xml")
    for host_guest_code in host_guest_codes:

        print(f"Preparing {host_guest_code}\n")
        try:
            system, model = prepare_openmm_system(host_guest_code, ff)
        except:
            print(f"ERROR trying to prepare {host_guest_code}\n")
            continue

        print(f"Simulating {host_guest_code}\n")
        try:
            simulate(host_guest_code, system, model)
        except:
            print(f"ERROR trying to simulate {host_guest_code}\n")
            continue
