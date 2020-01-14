#ifndef ZISA_MPI_PROGRESS_BAR_HPP_IUXXP
#define ZISA_MPI_PROGRESS_BAR_HPP_IUXXP

#include <zisa/config.hpp>
#include <zisa/io/progress_bar.hpp>
#include <zisa/parallelization/mpi.hpp>

namespace zisa {

class MPIProgressBar : public ProgressBar {
public:
  MPIProgressBar(std::shared_ptr<ProgressBar> serial_bar, MPI_Comm mpi_comm);

  void write_progress(std::ostream &out, const std::string &msg) override;

  void activate() override;
  void deactivate() override;
  void reset() override;

private:
  std::shared_ptr<ProgressBar> serial_bar;
  int rank;
};

}

#endif // ZISA_MPI_PROGRESS_BAR_HPP
