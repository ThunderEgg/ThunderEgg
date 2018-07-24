#ifndef GMGTriLinIntp_H
#define GMGTriLinIntp_H
#include "DomainCollection.h"
#include "InterLevelComm.h"
#include "Interpolator.h"
#include <memory>
namespace GMG
{
/**
 * @brief Simple class that directly places values from coarse vector to fine vector. (This is
 * O(0) accuracy, need to replace with bi-linear interpolation.
 */
class TriLinIntp : public Interpolator
{
	private:
	/**
	 * @brief The coarser set of domains
	 */
	std::shared_ptr<DomainCollection> coarse_dc;
	/**
	 * @brief The finer set of domains
	 */
	std::shared_ptr<DomainCollection> fine_dc;
	/**
	 * @brief The comm package between the levels.
	 */
	std::shared_ptr<InterLevelComm> ilc;

	public:
	/**
	 * @brief Create new TriLinIntp object.
	 *
	 * @param coarse_dc the coarser set of domains.
	 * @param fine_dc the finer set of domains.
	 * @param ilc the comm package between the levels.
	 */
	TriLinIntp(std::shared_ptr<DomainCollection> coarse_dc, std::shared_ptr<DomainCollection> fine_dc,
	         std::shared_ptr<InterLevelComm> ilc);
	/**
	 * @brief Interpolate from the finer level to the coarser level.
	 *
	 * @param coarse the input vector from the coarser level
	 * @param fine the output vector for the finer level
	 */
	void interpolate(PW<Vec> coarse, PW<Vec> fine) const;
};
} // namespace GMG
#endif
