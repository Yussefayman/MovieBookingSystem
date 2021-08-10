
package moviebookingsystem;

import database.mysql;
import java.io.IOException;
import java.net.URL;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.time.LocalDate;
import java.util.ResourceBundle;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Button;


public class Chair_bookingController implements Initializable {
   
    hall h = new hall();
    movie m= new movie();
    cinema c= new cinema();
    
    int mo,ci,ha,show_id1;
    @FXML public Button A1 ;
    @FXML public Button A2 ;
    @FXML public Button A3 ;
    @FXML public Button A4 ;
    @FXML public Button A5 ;
    @FXML public Button A6 ;
    
    @FXML public Button B1 ;
    @FXML public Button B2 ;
    @FXML public Button B3 ;
    @FXML public Button B4 ;
    @FXML public Button B5 ;
    @FXML public Button B6 ;
    
    @FXML public Button C1 ;
    @FXML public Button C2 ;
    @FXML public Button C3 ;
    @FXML public Button C4 ;
    @FXML public Button C5 ;
    @FXML public Button C6 ;
            
    @FXML
    public Button btn;
    String city= Submit_c_cController.city_;
    String cinema= Submit_c_cController.cinema_;
    String movie_name= FilmsController.movieName;
    String show_time = Date_time_submitController.time_;
    LocalDate show_date= Date_time_submitController.date_;
    String show_type = Date_time_submitController.type_;
    @FXML
    public void submit_chair(ActionEvent e) throws IOException{
        Controls control = new Controls();
         control.switchScenes("final_win");
    }
    @Override
    public void initialize(URL url, ResourceBundle rb) {
        try {
            check_chairs(movie_name,cinema ,city, show_time,show_type, show_date.toString());
        } catch (SQLException ex) {
            Logger.getLogger(Chair_bookingController.class.getName()).log(Level.SEVERE, null, ex);
        }
         
    }  
    @FXML
    public void submit_id(ActionEvent e) throws SQLException{
     
      btn = (Button) e.getSource();
        String id = btn.getId();
    if(h.hall_chairs(movie_name,cinema, city, id,show_time,show_type,show_date.toString() ))
          btn.setDisable(true);
      
    }
    @FXML
     public void  check_chairs(String movie_name_,String cinema_name,String cinema_city,String time,String show_type,String date) throws SQLException{
       mo=m.movie_id(movie_name_);
       ci=c.cinema_id(cinema_name, cinema_city);
       ha= h.hall_id(cinema_name,cinema_city);
       show_id1= h.show_id(time, show_type);
       ResultSet rs = mysql.RunOnQuery("SELECT * FROM chairscheck");
              while (rs.next()){
                 String chairid =rs.getString("chair_id");
                 int h_id =Integer.valueOf(rs.getString("hall_id"));
                 int movie = Integer.valueOf(rs.getString("movie_id"));
                 int cinema =Integer.valueOf(rs.getString("cinema_id"));
                 int s_id = Integer.valueOf(rs.getString("show_id"));
                 String type=rs.getString("show_type");
                 String s_date=rs.getString("date");
                 
                     if(chairid.equals("A1")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           A1.setDisable(true);
                          
                      }
                   else if(chairid.equals("A2")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           A2.setDisable(true);
                      }
                      else if(chairid.equals("A3")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           A3.setDisable(true);
                          
                      }
                      else if(chairid.equals("A4")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           A4.setDisable(true);
                          
                      }
                     else if(chairid.equals("A5")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           A5.setDisable(true);
                          
                      }
                     else if(chairid.equals("A6")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           A6.setDisable(true);
                          
                      }
                     else if(chairid.equals("B1")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           B1.setDisable(true);
                          
                      }
                    else if(chairid.equals("B2")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           B2.setDisable(true);
                            
                      }
                    else if(chairid.equals("B3")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           B3.setDisable(true);
                           
                      }
                    else if(chairid.equals("B4")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           B4.setDisable(true);
                           
                      }
                    else if(chairid.equals("B5")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           B5.setDisable(true);
                           
                      }
                    else if(chairid.equals("B6")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           B6.setDisable(true);
                            
                      }
                    else if(chairid.equals("C1")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           C1.setDisable(true);
                           
                      }
                    else if(chairid.equals("C2")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           C2.setDisable(true);
                            
                      }
                    else if(chairid.equals("C3")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           C3.setDisable(true);
                           
                      }
                    else if(chairid.equals("C4")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           C4.setDisable(true);
                           
                      }
                    else if(chairid.equals("C5")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           C5.setDisable(true);
                           
                      }
                    else if(chairid.equals("C6")&&h_id==ha&&movie==mo&&cinema==ci&&s_id==show_id1&&type.equals(show_type)&&s_date.equals(date)){
                           
                           C6.setDisable(true);
                            
                      }
                  }
                 
              }
   }

