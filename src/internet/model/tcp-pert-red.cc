#include "tcp-pert-red.h"
#include "ns3/log.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("TcpPertRed");
NS_OBJECT_ENSURE_REGISTERED (TcpPertRed);

TypeId
TcpPertRed::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::TcpPertRed")
    .SetParent<TcpNewReno> ()
    .AddConstructor<TcpPertRed> ()
    .SetGroupName ("Internet")
    .AddAttribute ("Thresh1", "Threshold 1",
                   TimeValue (Seconds(0.005)),
                   MakeTimeAccessor (&TcpPertRed::m_thresh1),
                  MakeTimeChecker ())
    .AddAttribute ("Thresh2", "Threshold 2",
                   TimeValue (Seconds(0.010)),
                   MakeTimeAccessor (&TcpPertRed::m_thresh2),
                  MakeTimeChecker ())
    .AddAttribute ("Thresh3", "Threshold 3",
                   TimeValue (Seconds(0.020)),
                   MakeTimeAccessor (&TcpPertRed::m_thresh3),
                  MakeTimeChecker ())
  ;
  return tid;
}

TcpPertRed::TcpPertRed (void)
  : TcpNewReno (),
    m_thresh1 (Time (Seconds (0.005))),
    m_thresh2 (Time (Seconds (0.010))),
    m_thresh3 (Time (Seconds (0.020)))
{
  NS_LOG_FUNCTION (this);
}

TcpPertRed::TcpPertRed (const TcpPertRed& sock)
  : TcpNewReno (sock),
    m_thresh1 (sock.m_thresh1),
    m_thresh2 (sock.m_thresh2),
    m_thresh3 (sock.m_thresh3)
	
{
  NS_LOG_FUNCTION (this);
}

TcpPertRed::~TcpPertRed (void)
{
  NS_LOG_FUNCTION (this);
}

Ptr<TcpCongestionOps>
TcpPertRed::Fork (void)
{
  return CopyObject<TcpPertRed> (this);
}

void
TcpPertRed::PktsAcked (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked,
                     const Time& rtt)
{
  NS_LOG_FUNCTION (this << tcb << segmentsAcked << rtt);

  if (rtt.IsZero ())
    {
      return;
    }

  m_minRtt = std::min (m_minRtt, rtt);
  NS_LOG_DEBUG ("Updated m_minRtt = " << m_minRtt);

  m_maxRtt = std::max (m_maxRtt, rtt);
  UpdatePertVars(rtt);
}

void
TcpPertRed::CongestionStateSet (Ptr<TcpSocketState> tcb,
                              const TcpSocketState::TcpCongState_t newState)
{
  NS_LOG_FUNCTION (this << tcb << newState);

}

void
TcpPertRed::IncreaseWindow (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked)
{
  NS_LOG_FUNCTION (this << tcb << segmentsAcked);
  NS_LOG_FUNCTION (this << tcb << segmentsAcked);
}

std::string
TcpPertRed::GetName () const
{
  return "TcpPertRed";
}

uint32_t
TcpPertRed::GetSsThresh (Ptr<const TcpSocketState> tcb,
                       uint32_t bytesInFlight)
{
  NS_LOG_FUNCTION (this << tcb << bytesInFlight);
  return std::max (uint32_t((1-m_beta)*tcb->m_ssThresh.Get ()), 2 * tcb->m_segmentSize); 
}
void TcpPertRed::UpdatePertVars(const Time& l_rtt)
{
}
void TcpPertRed::CheckChangeLossProb(Ptr<TcpSocketState> tcb)
{
}
void TcpPertRed::CalculateP ()
{
}
void TcpPertRed::CheckAndSetAlpha(Ptr<TcpSocketState> tcb)           
{
}
} // namespace ns3
