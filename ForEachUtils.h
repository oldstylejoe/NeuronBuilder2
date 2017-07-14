//Joe Snider
//6/06
//
//Utilities that take advantage of std::for_each.

#ifndef FOR_EACH_UTILS
#define FOR_EACH_UTILS

namespace feu {

	//Do averaging.
	//Call as:
	//     CAverage C1;
	//     for_each(begin, end, C1)
	//     mean = C1.GetAverage();
	//     standard_error = C1.GetStandardError();
	class CAverage{
	public:
		//initialize.
		CAverage(): m_dAverage(0.), m_dCount(0.) {}

        CAverage(const CAverage& inCopy) {
            *this = inCopy;
        }

        CAverage& operator=(const CAverage& inCopy) {
            m_dAverage = inCopy.m_dAverage;
            m_dCount = inCopy.m_dCount;
            m_dStdErr = inCopy.m_dStdErr;
            return *this;
        }

		//The function
		void operator()(const double& inX) {
			m_dStdErr += inX*inX;
			m_dAverage += inX;
			m_dCount += 1.;
		}

		//overload casting to get the average
		operator double() const {
			return GetAverage();
		}

		//Get the average.
		double GetAverage() const {
			if(m_dCount > 0.5) {
				return m_dAverage / m_dCount;
			}
			return 0.; //default return
		}

		//Get the std err.
		double GetStandardError() const {
			if(m_dCount > 0.5) {
				return sqrt(m_dStdErr - m_dAverage*m_dAverage/m_dCount) / m_dCount;
			}
			return 0.; //default return
		}

	private:
		double m_dAverage;
		double m_dCount;
		double m_dStdErr;
	};

}

#endif //FOR_EACH_UTILS
