/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */	
/*
 * Copyright (c) 2016 ResiliNets, ITTC, University of Kansas
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Truc Anh N. Nguyen <annguyen@ittc.ku.edu>
 *
 * James P.G. Sterbenz <jpgs@ittc.ku.edu>, director
 * ResiliNets Research Group  http://wiki.ittc.ku.edu/resilinets
 * Information and Telecommunication Technology Center (ITTC)
 * and Department of Electrical Engineering and Computer Science
 * The University of Kansas Lawrence, KS USA.
 */

#include "tcp-pert-pi.h"
#include "ns3/log.h"
#include "ns3/sequence-number.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("TcpPertPi");
NS_OBJECT_ENSURE_REGISTERED (TcpPertPi);

TypeId
TcpPertPi::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::TcpPertPi")
    .SetParent<TcpNewReno> ()
    .AddConstructor<TcpPertPi> ()
    .SetGroupName ("Internet")
    .AddAttribute ("Thresh1", "Threshold 1",
                   TimeValue (Seconds(0.005)),
                   MakeTimeAccessor (&TcpPertPi::m_thresh1),
                  MakeTimeChecker ())
    .AddAttribute ("Thresh2", "Threshold 2",
                   TimeValue (Seconds(0.010)),
                   MakeTimeAccessor (&TcpPertPi::m_thresh2),
                  MakeTimeChecker ())
    .AddAttribute ("Thresh3", "Threshold 3",
                   TimeValue (Seconds(0.020)),
                   MakeTimeAccessor (&TcpPertPi::m_thresh3),
                  MakeTimeChecker ())
    .AddAttribute ("MaxP", "Maximum Probability",
                   DoubleValue (0.050),
                   MakeDoubleAccessor (&TcpPertPi::m_maxp),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("Beta", "Pert Window Decrease Factor",
                   DoubleValue (0.350),
                   MakeDoubleAccessor (&TcpPertPi::m_beta),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("Alpha", "Pert Window Increase Factor",
                   DoubleValue (32),
                   MakeDoubleAccessor (&TcpPertPi::m_alpha),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("DropProbability", "Drop Probability",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertPi::m_dProb),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("EarlyResponseProbability", "Early Response Probability",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertPi::m_erProb),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("PertProbability", "Early Response Probability",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertPi::m_pertProb),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("ND", "Packets sent since last packet drop",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertPi::m_nd),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("NDP1", "Packets sent since last Early Response",
                   DoubleValue (0),
                   MakeDoubleAccessor (&TcpPertPi::m_ndp1),
                  MakeDoubleChecker <double> ())
    .AddAttribute ("QueueRef",
                   "Desired queue size",
                   DoubleValue (0.50),
                   MakeDoubleAccessor (&TcpPertPi::m_qRef),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("A",
                   "Value of A",
                   DoubleValue (0.00001822),
                   MakeDoubleAccessor (&TcpPertPi::m_a),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("B",
                   "Value of B",
                   DoubleValue (0.00001816),
                   MakeDoubleAccessor (&TcpPertPi::m_b),
                   MakeDoubleChecker<double> ())
   
  ;
  return tid;
}

TcpPertPi::TcpPertPi (void)
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
m_sender(false),
 m_qOld(0)
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
  m_rtrsEvent = Simulator::Schedule (Time (Seconds (0.3)), &TcpPertPi::CalculateP, this);
 
  NS_LOG_FUNCTION (this);
}

TcpPertPi::TcpPertPi (const TcpPertPi& sock)
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

TcpPertPi::~TcpPertPi (void)
{
  NS_LOG_FUNCTION (this);
}

Ptr<TcpCongestionOps>
TcpPertPi::Fork (void)
{
  return CopyObject<TcpPertPi> (this);
}

void
TcpPertPi::PktsAcked (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked,
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
  //printf("ACK : P : %f m_minRtt: %f t3 : %f\n", m_pertProb,m_minRtt.GetSeconds (),m_thresh3.GetSeconds());
  UpdatePertVars(rtt);
  //NS_LOG_DEBUG ("Updated m_cntRtt = " << m_cntRtt);
  if(!m_sender)
  {
        m_sender = true;
  }
}

void
TcpPertPi::CongestionStateSet (Ptr<TcpSocketState> tcb,
                              const TcpSocketState::TcpCongState_t newState)
{
  NS_LOG_FUNCTION (this << tcb << newState);
        std::cout<<"statechange:::::::::::::::::::::::\n";
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
  printf("CAR : cwnd : %d sst : %d \n", tcb->m_cWnd.Get(),tcb->m_ssThresh.Get());
  }

}

void
TcpPertPi::IncreaseWindow (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked)
{
  NS_LOG_FUNCTION (this << tcb << segmentsAcked);

  NS_LOG_FUNCTION (this << tcb << segmentsAcked);
  bool early_response = false;
  Ptr<UniformRandomVariable> prob_this_pkt = CreateObject<UniformRandomVariable> ();  
//proactive_slowdown 
        
 //printf("IW : P : %f m_cWnd: %d m_sst: %d prob_this_pkt : %f thisS: %d lastS : %d\n", m_pertProb,tcb->m_cWnd.Get (),tcb->m_ssThresh.Get(),prob_this_pkt->GetValue() ,tcb->m_lastAckedSeq.GetValue(),m_lastEarlyResponseSeq.GetValue());
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
  
    tcb->m_cWnd = (1-m_beta) * tcb->m_cWnd.Get ();
    tcb->m_ssThresh = GetSsThresh (tcb, 0);
    early_response = true;
    m_lastEarlyResponseSeq = tcb->m_lastAckedSeq;
       printf("PS : cwnd : %d sst : %d \n", tcb->m_cWnd.Get(),tcb->m_ssThresh.Get());
        std::cout<<m_lastEarlyResponseSeq;
    }


        
  if (!early_response)
    {
      m_nd = m_nd + 1;
      m_ndp1 = m_ndp1 + 1;
      CheckChangeLossProb(tcb);
      if (tcb->m_cWnd < tcb->m_ssThresh)
        {
           TcpNewReno::SlowStart (tcb, segmentsAcked);
           //printf("SS : cwnd : %d sst : %d \n", tcb->m_cWnd.Get(),tcb->m_ssThresh.Get());
        }

      if (tcb->m_cWnd >= tcb->m_ssThresh)
        {
          double adder = static_cast<double> (m_alpha) / tcb->m_cWnd.Get ();
          tcb->m_cWnd += static_cast<uint32_t> (adder);
          NS_LOG_INFO ("In CongAvoid, updated to cwnd " << tcb->m_cWnd <<
                       " ssthresh " << tcb->m_ssThresh);
          //printf("CA : cwnd : %d sst : %d \n", tcb->m_cWnd.Get(),tcb->m_ssThresh.Get());
        }
    }
  

}

std::string
TcpPertPi::GetName () const
{
  return "TcpPertPi";
}

uint32_t
TcpPertPi::GetSsThresh (Ptr<const TcpSocketState> tcb,
                       uint32_t bytesInFlight)
{
  NS_LOG_FUNCTION (this << tcb << bytesInFlight);
        //printf("SST : sst : %d\n",std::max (static_cast<uint32_t>((1-m_beta)*tcb->m_cWnd.Get ()), 2 * tcb->m_segmentSize));
  return std::max (static_cast<uint32_t>((1-m_beta)*tcb->m_cWnd.Get ()), 2 * tcb->m_segmentSize); //need to modify according to pert
}









//printf("IW : m_cWnd : %d",std::max (static_cast<uint32_t>((1-m_beta)*tcb->m_cWnd.Get ()), 2 * tcb->m_segmentSize));


void TcpPertPi::UpdatePertVars(const Time& l_rtt)
{
	m_thresh3 = Time(Seconds (std::max ((2*m_thresh2.GetSeconds ()), 0.65*(m_maxRtt.GetSeconds () - m_minRtt.GetSeconds ()))));
	if (m_pertSrtt> 0) {
		m_pertSrtt = m_pertSrtt*0.99 + 0.01 *l_rtt.GetSeconds ();
	}
	if (m_pertSrtt <= 0) {
	        m_pertSrtt = l_rtt.GetSeconds ();            
	}      
	double maxq = std::max(0.010, (m_maxRtt.GetSeconds () - m_minRtt.GetSeconds ())); ////Type Conversion //done
	double curq = m_pertSrtt - m_minRtt.GetSeconds ();
	if (curq > 0)
		m_beta = (curq) / ((curq + maxq));  ////Type conversion
          //printf("UPV P : %f m_pertSrtt: %f m_minRtt: %f beta : %f\n", m_pertProb,m_pertSrtt,m_minRtt.GetSeconds(),m_beta);
}


















void TcpPertPi::CheckChangeLossProb(Ptr<TcpSocketState> tcb)
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
          //printf("CLP : P : %f m_minRtt: %f m_dProb : %f\n", m_pertProb,m_minRtt.GetSeconds(),m_dProb);
}










void TcpPertPi::CalculateP ()
{
         NS_LOG_FUNCTION (this);
        double curr_srtt = m_pertSrtt;
        double p = 0.0;
        double minRtt = m_minRtt.GetSeconds ();
	double curq = curr_srtt - minRtt;
        
        {
          p = m_a * (curq - m_qRef) - m_b * (m_qOld - m_qRef) + m_pertProb;
        }
        p = (p < 0) ? 0 : p;
        p = (p > 1) ? 1 : p;

        m_qOld = curq;
	m_pertProb = p;
        m_rtrsEvent = Simulator::Schedule (Time (Seconds (1.0/170.0)), &TcpPertPi::CalculateP, this);
    
}



void TcpPertPi::CheckAndSetAlpha(Ptr<TcpSocketState> tcb)           ////********changed function name.
{

	double curr_rtt = m_pertSrtt;              
	double maxq = std::max(0.010, (m_maxRtt.GetSeconds () - m_minRtt.GetSeconds ())); ////Type Conversion
	double curq = curr_rtt - m_minRtt.GetSeconds ();
	double k1,k2,pp1,target;
	//first determine the region we are operating in: high speed, safe, or competitive
	if (curq <= m_thresh1.GetSeconds ())
	{
		//then we are in the high speed region
		m_highspeedRegionCounter++;               ////******add m_highSpeedRegionCounter in header file
		m_mode=0;                                  //// ****add m_mode in header file
	}
	else if (curq >= 0.5*maxq)
	{
		//then we are in the competitive region
		m_competeRegionCounter++;					////*******add m_competeRegionCounter in header file
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
        //// can we replace Scheduler::instance().clock() by simulator::now()


		//if we are in the compete region then alpha=
		if (m_competeRegionCounter >= tcb->m_cWnd)
		{
			//if the drop probability is not zero, then
			if (m_dProb != 0)
			{
				pp1=(1 + (m_erProb / m_dProb));
				if (curq > (maxq/2) && curq < (0.65*maxq)) {
					k1 = ((pp1-1)*maxq*16)/(15*100);
					k2 = 1 + ((k1*100)/maxq);
					target = k2 - (k1/(curq - 0.49*maxq));
				}
				else
					target = pp1;

				m_alpha = std::min( m_alpha+0.1, target);
				m_alpha = std::min(m_alpha, m_alphaMax);
//				//printf("\n Flow: %d, compete region target: %lf, max: %lf\n", flow_id, target, pp1);
			}
			//else, if the drop probability is zero:
			else
			{
				//verify: replaced 10*PRECISION with 10 in the equation -- according to discussion with Kiran
				//we can remove PRECISION, as we have floating point support in ns2 -- Kiran Kotla
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
//			printf("\n Flow: %d, safe region target: %lf\n", flow_id, target);
			//else
		}
		//now, reset the counters
		m_competeRegionCounter=0;
		m_highspeedRegionCounter=0;
		//update the last time alpha was updated
		m_lastUpdatedAlpha=Simulator::Now ();  //// we have to replace this line ////done
	//	if (m_flowId == 1)
		//printf("\n %lf: %lf %lf %lf %lf %lf %lf", Scheduler::instance().clock(), minRtt, curr_rtt, maxRtt, m_erProb, m_dProb, m_alpha);

	}
          //printf("CSA : P : %f curq: %f m_minRtt: %f m_alpha : %f\n", m_pertProb,curq,m_minRtt.GetSeconds(),m_alpha);
}







} // namespace ns3
