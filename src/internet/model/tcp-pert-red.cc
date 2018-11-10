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
                  MakeTimeChecker ()),
		.AddAttribute ("MaxP", "Maximum Probability",
                   DoubleValue (0.050),
                   MakeDoubleAccessor (&TcpPertRed::m_maxp),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("Beta", "Pert Window Decrease Factor",
                   DoubleValue (0.350),
                   MakeDoubleAccessor (&TcpPertRed::m_beta),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("Alpha", "Pert Window Increase Factor",
                   DoubleValue (32),
                   MakeDoubleAccessor (&TcpPertRed::m_alpha),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("DropProbability", "Drop Probability",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertRed::m_dProb),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("EarlyResponseProbability", "Early Response Probability",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertRed::m_erProb),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("PertProbability", "Early Response Probability",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertRed::m_pertProb),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("ND", "Packets sent since last packet drop",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertRed::m_nd),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("NDP1", "Packets sent since last Early Response",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertRed::m_ndp1),
                  MakeDoubleChecker <double> ())
  ;
  return tid;
}

TcpPertRed::TcpPertRed (void)
  : TcpNewReno (),
    m_thresh1 (Time (Seconds (0.005))),
    m_thresh2 (Time (Seconds (0.010))),
    m_thresh3 (Time (Seconds (0.020))),
    m_maxp (0.050),
    m_beta (0.35),       
    m_alpha (33), 
    m_dProb (0),
    m_erProb (0),
    m_pertProb (0),
    m_nd (0),
    m_ndp1 (0),
    m_minRtt(Time(Seconds(0.05))),
 	 	m_maxRtt(Time(Seconds(0.00))),
		m_pertSrtt(0),
  	m_changeWindow(0),
  	m_competeRegionCounter(0),
  	m_highspeedRegionCounter(0),
  	m_mode(0),
  	m_maxProb(1),
m_lastUpdatedAlpha(0)
{
  m_weight.push_back(0.2);
  m_weight.push_back(0.4);
  m_weight.push_back(0.6);
  m_weight.push_back(0.8);
  m_weight.push_back(1);
  m_weight.push_back(1);
  m_weight.push_back(1);
  m_weight.push_back(1);
  for(int i = 0 ; i < 8 ; i++)
  {
    m_historyND.push_back(0);
    m_historyNDp1.push_back(0);
  }
  m_rtrsEvent = Simulator::Schedule (Time (Seconds (1.0 / 170)), &TcpPertRed::CalculateP, this);
  NS_LOG_FUNCTION (this);
}

TcpPertRed::TcpPertRed (const TcpPertRed& sock)
  : TcpNewReno (sock),
    m_thresh1 (sock.m_thresh1),
    m_thresh2 (sock.m_thresh2),
    m_thresh3 (sock.m_thresh3),
    m_maxp (sock.m_maxp),
    m_beta (sock.m_beta),
    m_alpha (sock.m_alpha ),
    m_dProb (sock.m_dProb),
    m_erProb (sock.m_erProb),
    m_pertProb (sock.m_pertProb),
    m_nd (sock.m_nd),
    m_ndp1 (sock.m_ndp1),
    m_minRtt (sock.m_minRtt),
    m_maxRtt (sock.m_maxRtt),
    m_pertSrtt(sock.m_pertSrtt),
    m_changeWindow(sock.m_changeWindow),
		m_competeRegionCounter(sock.m_competeRegionCounter),
		m_highspeedRegionCounter(sock.m_highspeedRegionCounter),
		m_mode(sock.m_mode),  
  	m_maxProb(1),
		m_lastUpdatedAlpha(sock.m_lastUpdatedAlpha)
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
  m_maxRtt = std::max (m_maxRtt, rtt);
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
	bool early_response = false;
  Ptr<UniformRandomVariable> prob_this_pkt = CreateObject<UniformRandomVariable> ();  
	//proactive_slowdown 
  if (tcb->m_cWnd > tcb->m_ssThresh && (prob_this_pkt->GetValue() <= m_pertProb) && (tcb->m_lastAckedSeq >= m_lastEarlyResponseSeq))
  {
		if (m_ndp1 != 0)
    {
			m_historyNDp1.erase(m_historyNDp1.begin());
      m_historyNDp1.push_back(m_ndp1);
		  double den = 0;
		  for (int i=0;i<8;i++)
      {
			  den = den + m_weight[i]*m_historyNDp1[i];
      }		    
      m_erProb = 6/(den);
		  m_ndp1 = 0;                                         
	  }
 		tcb->m_cWnd = (1-m_beta) * tcb->m_cWnd.Get ();
    tcb->m_ssThresh = GetSsThresh (tcb, 0);
    early_response = true;
    m_lastEarlyResponseSeq = tcb->m_lastAckedSeq;
  }      
  if (!early_response)
  {
    m_nd = m_nd + 1;
    m_ndp1 = m_ndp1 + 1;
    CheckChangeLossProb(tcb);
    if (tcb->m_cWnd < tcb->m_ssThresh)
    {
    	TcpNewReno::SlowStart (tcb, segmentsAcked);
    }
		if (tcb->m_cWnd >= tcb->m_ssThresh)
    {
    	double adder = static_cast<double> (m_alpha) / tcb->m_cWnd.Get ();
      tcb->m_cWnd += static_cast<uint32_t> (adder);
    }
	}
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
