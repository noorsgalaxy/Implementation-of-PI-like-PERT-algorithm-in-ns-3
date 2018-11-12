#include "tcp-pert-red.h"
#include "ns3/log.h"
#include "ns3/sequence-number.h"

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
    m_minRtt(Time::Max ()),
  m_maxRtt(Time::Min ()),
  m_pertSrtt(0),
  m_changeWindow(0),
  m_competeRegionCounter(0),
  m_highspeedRegionCounter(0),
  m_mode(0),
  m_maxProb(1),
  m_lastUpdatedAlpha(0),
  m_alphaMax(32),
  m_sender(false)
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
  m_rtrsEvent = Simulator::Schedule (Time (Seconds (0.3)), &TcpPertRed::CalculateP, this);
 
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
    m_changeWindow(sock.m_changeWindow),
	m_competeRegionCounter(sock.m_competeRegionCounter),
	m_highspeedRegionCounter(sock.m_highspeedRegionCounter),
	m_mode(sock.m_mode),  
  m_maxProb(sock.m_maxProb),
	m_lastUpdatedAlpha(sock.m_lastUpdatedAlpha),
m_alphaMax(sock.m_alphaMax)
	
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
  if(newState ==  TcpSocketState::CA_RECOVERY or newState ==  TcpSocketState::CA_LOSS )
  {
    m_changeWindow = 0;
	  if (m_nd != 0) {
		  m_historyND.erase(m_historyND.begin());
		  m_historyND.push_back(m_nd);
		  double L = 0;
		  for (int i=0;i<8;i++)
      {
		    L = L + m_weight[i]*m_historyND[i];
      }		  
      m_dProb = 6/(L);
		  m_nd = 0; 
	}
   }
}

void
TcpPertRed::IncreaseWindow (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked)
{
  NS_LOG_FUNCTION (this << tcb << segmentsAcked);

  NS_LOG_FUNCTION (this << tcb << segmentsAcked);
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
		    m_ndp1 = 0;                                         //Reset the epoch packets
	    }
    tcb->m_cWnd = ((1-m_beta) * tcb->GetCwndInSegments ()) * tcb->m_segmentSize;
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
          double adder = static_cast<double> (m_alpha) / tcb->GetCwndInSegments ();
          tcb->m_cWnd += static_cast<uint32_t> (adder)* tcb->m_segmentSize;
          NS_LOG_INFO ("In CongAvoid, updated to cwnd " << tcb->m_cWnd <<
                       " ssthresh " << tcb->m_ssThresh);
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
       
  return std::max (static_cast<uint32_t>(((1-m_beta) * tcb->GetCwndInSegments ()) * tcb->m_segmentSize), 2 * tcb->m_segmentSize);
}

void TcpPertRed::UpdatePertVars(const Time& l_rtt)
{
	m_thresh3 = Time(Seconds (std::max ((2*m_thresh2.GetSeconds ()), 0.65*(m_maxRtt.GetSeconds () - m_minRtt.GetSeconds ()))));
	if (m_pertSrtt> 0) 
  {
		m_pertSrtt = m_pertSrtt*0.99 + 0.01 *l_rtt.GetSeconds ();
	}
	if (m_pertSrtt <= 0) 
  {
	        m_pertSrtt = l_rtt.GetSeconds ();            
	}      
	double maxq = std::max(0.010, (m_maxRtt.GetSeconds () - m_minRtt.GetSeconds ())); 
	double curq = m_pertSrtt - m_minRtt.GetSeconds ();
	if (curq > 0)
		m_beta = (curq) / ((curq + maxq)); 
    }

void TcpPertRed::CheckChangeLossProb(Ptr<TcpSocketState> tcb)
{
	double L = 0;
	for (int i=0;i<8;i++)
		L = L + m_weight[i]*m_historyND[i];
	L = L/6;
	if (m_nd < L)
		return;
	if (m_changeWindow == 0)
	{
		m_changeWindow = 1;                   
		m_historyND.erase(m_historyND.begin());
	}
	m_historyND.push_back(m_nd);
	L = 0;
	for (int i=0;i<8;i++)
		L = L + m_weight[i]*m_historyND[i];
	m_dProb = 6/(L);
	CheckAndSetAlpha (tcb);     
}

void TcpPertRed::CalculateP ()
{
	double curr_srtt = m_pertSrtt;
	double p = m_pertProb;
        double minRtt = m_minRtt.GetSeconds ();
	double curq = curr_srtt - minRtt;
	m_thresh3 = Time(Seconds (std::max ((2*m_thresh2.GetSeconds ()), 0.65*(m_maxRtt.GetSeconds () - m_minRtt.GetSeconds ()))));

	if (curr_srtt >= (minRtt + m_thresh2.GetSeconds ()))
	{
		p = m_maxp + ((m_maxProb-m_maxp)*((curq - m_thresh2.GetSeconds ())) / (m_thresh3.GetSeconds () - m_thresh2.GetSeconds ()));
	}
	else if (curr_srtt >= (minRtt + m_thresh1.GetSeconds ()))
	{
      		p = m_maxp *((curq - m_thresh1.GetSeconds ())/(m_thresh2.GetSeconds () - m_thresh1.GetSeconds ()));
	}
if(m_sender)
  printf("P : %f curr_srtt: %f curq: %f m_minRtt: %f t3 : %f\n", m_pertProb,curr_srtt,curq,m_minRtt.GetSeconds(),m_thresh3.GetSeconds());
	if (p < 0) p = 0.0;
	if (p > 1) p = 1;
	m_pertProb = p;
    m_rtrsEvent = Simulator::Schedule (Time (Seconds (1.0/170.0)), &TcpPertRed::CalculateP, this);
}

void TcpPertRed::CheckAndSetAlpha(Ptr<TcpSocketState> tcb)           ////********changed function name.
{
	double curr_rtt = m_pertSrtt;              
	double maxq = std::max(0.010, (m_maxRtt.GetSeconds () - m_minRtt.GetSeconds ())); ////Type Conversion
	double curq = curr_rtt - m_minRtt.GetSeconds ();
	double k1,k2,pp1,target;
	if (curq <= m_thresh1.GetSeconds ())
	{
		//then we are in the high speed region
		m_highspeedRegionCounter++;               
		m_mode=0;                                  
	}
	else if (curq >= 0.5*maxq)
	{
		//then we are in the competitive region
		m_competeRegionCounter++;					
		m_mode = 2;
	}
	else
	{
		//else, we are in the safe region
		m_highspeedRegionCounter = 0;
		m_competeRegionCounter = 0;
		m_mode = 1;
	}
	//now, update alpha accordingly
	//if 5 rtts have elapsed, then update alpha accordingly
	if (Simulator::Now ().GetSeconds () - m_lastUpdatedAlpha.GetSeconds () >= 5*curr_rtt)        /////Time Type
	{
		//if we are in the compete region then alpha=
		if (m_competeRegionCounter >= tcb->m_cWnd)
		{
			//if the drop probability is not zero, then
			if (m_dProb != 0)
			{
				pp1=(1 + (m_erProb / m_dProb));
				if (curq > (maxq/2) && curq < (0.65*maxq)) 
        {
					k1 = ((pp1-1)*maxq*16)/(15*100);
					k2 = 1 + ((k1*100)/maxq);
					target = k2 - (k1/(curq - 0.49*maxq));
				}
				else
					target = pp1;

				m_alpha = std::min( m_alpha+0.1, target);
				m_alpha = std::min(m_alpha, m_alphaMax);
			}
			//else, if the drop probability is zero:
			else
			{
      		m_alpha = std::min(m_alpha+0.1, m_alphaMax);
			}
		}
		// if we are in the highspeed region then alpha=
		else if	(m_highspeedRegionCounter >= tcb->m_cWnd)
		{
			//update alpha
			m_alpha = std::min((m_alpha+0.5), m_alphaMax);
		}
		//else, if we are in the safe region, then alpha=
		else
		{
			//if the current queue is between the min threshold and the max threshold
			if (curq > m_thresh1.GetSeconds () && curq < m_thresh2.GetSeconds ()) {
				k2 = (m_thresh2.GetSeconds () - m_thresh1.GetSeconds ())/31;
				k1 = k2 + m_thresh2.GetSeconds ();
				target = (k1)/(k2+curq);
			}
			else
				target = 1;

			m_alpha = 0.9*m_alpha;

			if (m_alpha < 1)
				m_alpha = 1;
		}
		//now, reset the counters
		m_competeRegionCounter=0;
		m_highspeedRegionCounter=0;
		//update the last time alpha was updated
		m_lastUpdatedAlpha=Simulator::Now (); 
	 }
  }
} // namespace ns3
