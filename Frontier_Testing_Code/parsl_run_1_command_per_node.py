'''
parsl_run_1_command_per_node.py

The purpose of this script is to create a way for executing a file of commands
that each need to run on a single node on Frontier.

***This needs to be run in the parsl_bash conda environment ***
'''

import argparse
import parsl
from parsl.app.app import bash_app
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.addresses import address_by_hostname
import os

# Get the number of nodes:
node_raw = os.getenv("SLURM_NODELIST")
node_list = [f"frontier{x}" for x  in node_raw.lstrip("frontier[").rstrip("]").split(",")]
num_nodes = len(node_list)

# Parsl configuration for OLCF Frontier (56 cores per node)
config = Config(
    executors=[
        HighThroughputExecutor(
            label="frontier_htex",
            address=address_by_hostname(),
            cores_per_worker=1.0,  # One worker per task
            max_workers_per_node=1,  # One worker per node
            provider=LocalProvider(
                # Number of nodes job
                nodes_per_block=num_nodes,
                launcher=SrunLauncher(overrides='-c 56'),
                init_blocks=1,
                max_blocks=1,
            ),
        )
    ],
    usage_tracking=True,
)

parsl.clear()
parsl.load(config)

@bash_app
def run_command(cmd):
    return f"""
    {cmd}
    """

def main():
    parser = argparse.ArgumentParser(
        description='Run unix commands in parallel on Frontier using Parsl'
    )
    parser.add_argument(
        'command_file',
        type=str,
        help='Path to text file with unix commands, one per line'
    )
    args = parser.parse_args()

    # Read commands from file
    with open(args.command_file, 'r') as f:
        commands = [line.strip() for line in f if line.strip()]

    # Submit each command as a Parsl task
    tasks = []
    for cmd in commands:
        task = run_command(cmd)
        tasks.append(task)

    # Wait for all tasks to complete
    results = []
    for i, task in enumerate(tasks):
        try:
            result = task.result()
            print(f"SUCCESS: Task {i} completed successfully")
            results.append(result)
        except Exception as e:
            print(f"ERROR: Task {i} failed with error: {e}")
            results.append(None)

    parsl.clear()

if __name__ == '__main__':
    main()
