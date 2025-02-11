#include "../../../Block_Cost.hxx"
#include "../../../Block_Map.hxx"
#include "../../Block_Info.hxx"

std::vector<std::vector<Block_Map>>
compute_block_grid_mapping(const size_t &procs_per_node,
                           const size_t &num_nodes,
                           const std::vector<Block_Cost> &block_costs);

void Block_Info::allocate_blocks(const std::vector<Block_Cost> &block_costs,
                                 const size_t &procs_per_node,
                                 const size_t &proc_granularity,
                                 const Verbosity &verbosity)
{
  // Reverse sort, with largest first
  std::vector<Block_Cost> sorted_costs(block_costs);
  std::sort(sorted_costs.rbegin(), sorted_costs.rend());
  const size_t num_procs(El::mpi::Size(El::mpi::COMM_WORLD));
  if(num_procs % procs_per_node != 0)
    {
      throw std::runtime_error(
        "Incompatible number of MPI processes and processes per node.  "
        "procsPerNode must evenly divide to total number of MPI "
        "processes:\n\tMPI processes: "
        + std::to_string(num_procs)
        + "\n\tprocsPerNode: " + std::to_string(procs_per_node));
    }
  if(procs_per_node % proc_granularity != 0)
    {
      throw std::runtime_error(
        "Incompatible number of processes per node and process granularity.  "
        "procGranularity mush evenly divide procsPerNode:\n\tprocsPerNode: "
        + std::to_string(procs_per_node)
        + "\n\tprocGranularity: " + std::to_string(proc_granularity));
    }
  const size_t num_nodes(num_procs / procs_per_node);
  std::vector<std::vector<Block_Map>> mapping(compute_block_grid_mapping(
    procs_per_node / proc_granularity, num_nodes, sorted_costs));

  for(auto &block_vector : mapping)
    for(auto &block_map : block_vector)
      {
        block_map.num_procs *= proc_granularity;
      }

  // Create an mpi::Group for each set of processors.
  El::mpi::Group default_mpi_group;
  El::mpi::CommGroup(El::mpi::COMM_WORLD, default_mpi_group);

  int rank(El::mpi::Rank(El::mpi::COMM_WORLD));

  if(verbosity >= Verbosity::debug && rank == 0)
    {
      std::stringstream ss;
      ss << "Block Grid Mapping\n"
         << "Node\tNum Procs\tCost\t\tBlock List\n"
         << "==================================================\n";
      for(size_t node = 0; node < mapping.size(); ++node)
        {
          for(auto &m : mapping[node])
            {
              ss << node << "\t" << m.num_procs << "\t\t"
                 << m.cost / static_cast<double>(m.num_procs) << "\t{";
              for(size_t ii = 0; ii < m.block_indices.size(); ++ii)
                {
                  if(ii != 0)
                    {
                      ss << ",";
                    }
                  ss << "(" << m.block_indices[ii] << ","
                     << num_points[m.block_indices[ii]] << ")";
                }
              ss << "}\n";
            }
          ss << "\n";
        }
      El::Output(ss.str());
    }

  int rank_begin(0), rank_end(0);
  for(auto &block_vector : mapping)
    {
      for(auto &block_map : block_vector)
        {
          rank_begin = rank_end;
          rank_end += block_map.num_procs;
          if(rank_end > rank)
            {
              block_indices = block_map.block_indices;
              break;
            }
        }
      if(rank_end > rank)
        {
          break;
        }
    }
  // We should be generating blocks to cover all of the processors,
  // even if there are more nodes than procs.  So this is a sanity
  // check in case we messed up something in
  // compute_block_grid_mapping.
  if(rank_end <= rank)
    {
      throw std::runtime_error("INTERNAL ERROR: Some procs were not covered "
                               "by compute_block_grid_mapping.\n"
                               "\trank = "
                               + std::to_string(rank) + "\n"
                               + "\trank_end = " + std::to_string(rank_end));
    }
  {
    std::vector<int> group_ranks(rank_end - rank_begin);
    std::iota(group_ranks.begin(), group_ranks.end(), rank_begin);
    El::mpi::Incl(default_mpi_group, group_ranks.size(), group_ranks.data(),
                  mpi_group.value);
  }
  El::mpi::Create(El::mpi::COMM_WORLD, mpi_group.value, mpi_comm.value);
}
