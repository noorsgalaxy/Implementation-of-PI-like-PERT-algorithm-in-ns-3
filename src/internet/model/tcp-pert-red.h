#ifndef TCPPERTRED_H
#define TCPPERTRED_H

#include "ns3/nstime.h"
#include "ns3/boolean.h"
#include "ns3/data-rate.h"
#include "ns3/timer.h"
#include "ns3/event-id.h"
#include "ns3/random-variable-stream.h"
#include "ns3/tcp-congestion-ops.h"
#include "ns3/tcp-recovery-ops.h"

namespace ns3 {

class TcpPertRed : public TcpNewReno
{
public:

  static TypeId GetTypeId (void);


  TcpPertRed (void);


  TcpPertRed (const TcpPertRed& sock);
  virtual ~TcpPertRed (void);

  virtual std::string GetName () const;


  virtual void PktsAcked (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked,
                          const Time& rtt);


  virtual void CongestionStateSet (Ptr<TcpSocketState> tcb,
                                   const TcpSocketState::TcpCongState_t newState);

  virtual void IncreaseWindow (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked);


  virtual uint32_t GetSsThresh (Ptr<const TcpSocketState> tcb,
                                uint32_t bytesInFlight);

  virtual Ptr<TcpCongestionOps> Fork ();


private:

  void CheckChangeLossProb(Ptr<TcpSocketState> tcb);
  void CalculateP();
  void CheckAndSetAlpha(Ptr<TcpSocketState> tcb);
  void UpdatePertVars(const Time& l_rtt);
  Time m_thresh1;                   //!<
  Time m_thresh2;                   //!<
  Time m_thresh3;                   //!<
  double m_maxp;                  //!<probability value at thresh2_
  double m_beta;                  //!< Decresing Factor
  double m_alpha;                 //!<Maximum alpha value set to 32
  double m_dProb;                 //!<weighted average of drop probability (p) computed in TFRC style
  double m_erProb;                //!<weighted average of early response probability (p') computed in TFRC style
  double m_pertProb;                //!<
  double m_nd;                    //!<number of packets sent in the current since last packet loss
  double m_ndp1;                  //!<number of packets sent since last early response (for the estimation of p')
  Time m_minRtt;  
  Time m_maxRtt;
  double m_pertSrtt;  
  uint32_t m_changeWindow;
  uint32_t m_competeRegionCounter;
  uint32_t  m_highspeedRegionCounter;
  uint32_t  m_mode;
  double m_maxProb;   
  Time  m_lastUpdatedAlpha;     
  double m_alphaMax;
  bool m_sender;
  double m_this;
  EventId m_rtrsEvent;
  SequenceNumber32 m_lastEarlyResponseSeq;
  std::vector<double> m_weight;
  std::vector<double> m_historyND;
  std::vector<double> m_historyNDp1;
};

} // namespace ns3

#endif // TCPPERTRED_H
