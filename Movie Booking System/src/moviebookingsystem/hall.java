package moviebookingsystem;

import database.mysql;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import javafx.scene.control.Alert;
import javax.swing.JOptionPane;


public class hall {
   int id;
   public static ArrayList<String> chairs = new ArrayList<String>();
   public hall(int id){
       this.id=id;
   }
   public hall(){
       
   }
    public int hall_id(String cinema_name,String cinema_city) throws SQLException{
         cinema c= new cinema();
        int x= c.cinema_id(cinema_name, cinema_city);
            ResultSet rs = mysql.RunOnQuery("SELECT * FROM hall");
              while (rs.next()){
                   id = Integer.valueOf(rs.getString("hall_id"));
                  int cinema = Integer.valueOf(rs.getString("cinema_id"));
                 if(cinema==x){
                     return id;
                  }
              }
        return 0;
           
        
   }
   public boolean check_hall(int hall_id,int cinema_id) throws SQLException{
       boolean found=false;
        ResultSet rs= mysql.RunOnQuery("select * from hall");
       while(rs.next()){
           int hall_id_= Integer.valueOf(rs.getString("hall_id"));
           int cin_id= Integer.valueOf(rs.getString("cinema_id"));
          if(hall_id_==hall_id&&cin_id==cinema_id){
              found=true;
           }
       }
       if(found==true){
           Alert nameAlert = new Alert(Alert.AlertType.WARNING);
                nameAlert.setTitle("WARNING: Adding Failed");
                nameAlert.setContentText("Hall and cinema are already exist :( ");
                nameAlert.showAndWait();
                return false; 
       }
       else{
     
         return true;
   }
       
}
   public void save_hall(int hall_id,int cinema_id,String movie_name) throws SQLException{
       movie m= new movie();
       int m_id= m.movie_id(movie_name);
       mysql.runnonquery("insert into hall (hall_id,cinema_id,movie_id) values ('"+hall_id+"','"+cinema_id+"','"+m_id+"')");
         
   }
 public void save_show(String cinema_name,String cinema_city,String movie_name) throws SQLException{
       cinema c= new cinema();
       movie m =new movie();
       String t1="2D";
       String t2="3D";
       String t3="4D";
      int h_id= hall_id(cinema_name,cinema_city);
       int c_id=c.cinema_id(cinema_name, cinema_city);
       int m_id= m.movie_id(movie_name);
       for(int i=1;i<=5;i++){
         if(i==1){
              mysql.runnonquery("insert into hall_has_show (hall_id,cinema_id,show_id,show_type,movie_id) values ('"+h_id+"','"+c_id+"','"+i+"','"+t1+"','"+m_id+"')");
            
         }  
         else if(i==2){
             mysql.runnonquery("insert into hall_has_show (hall_id,cinema_id,show_id,show_type,movie_id) values ('"+h_id+"','"+c_id+"','"+i+"','"+t1+"','"+m_id+"')");
             
         }
         else if(i==3){
            mysql.runnonquery("insert into hall_has_show (hall_id,cinema_id,show_id,show_type,movie_id) values ('"+h_id+"','"+c_id+"','"+i+"','"+t1+"','"+m_id+"')");
              
         }
         else if(i==4){
              mysql.runnonquery("insert into hall_has_show (hall_id,cinema_id,show_id,show_type,movie_id) values ('"+h_id+"','"+c_id+"','"+i+"','"+t2+"','"+m_id+"')");
         }
         else if(i==5){
              mysql.runnonquery("insert into hall_has_show (hall_id,cinema_id,show_id,show_type,movie_id) values ('"+h_id+"','"+c_id+"','"+i+"','"+t3+"','"+m_id+"')");
            
         }
       }
   }
 public boolean  hall_chairs(String movie_name,String cinema_name,String cinema_city,String click,String time,String show_type,String date) throws SQLException{
         
          movie m= new movie();
          cinema c = new cinema();
       int movie=   m.movie_id(movie_name);
       int cinema= c.cinema_id(cinema_name, cinema_city);
       int hall= hall_id(cinema_name,cinema_city);
       int show_id= show_id(time, show_type);
        mysql.runnonquery("insert into chairscheck (chair_id,hall_id,movie_id,cinema_id,show_id,show_type,date) values ('"+click+"','"+hall+"','"+movie+"','"+cinema+"','"+show_id+"','"+show_type+"','"+date+"')");
         chairs.add(click);
        
        Alert pAlert = new Alert(Alert.AlertType.INFORMATION);
                           pAlert.setTitle("Confirmation");
                           pAlert.setContentText("You Booked This Chair Successfully");
                           pAlert.showAndWait();
        return true;   
         }

 public int show_id(String time,String show_type) throws SQLException{
     
       if(time.equals("15:00:00") && show_type.equals("2D"))
       {
           
           return 1;
       }
       else if(time.equals("18:00:00") && show_type.equals("2D")){
           return 2;
       }
       else if(time.equals("21:00:00") && show_type.equals("2D")){
           return 3;
       }
       else if(time.equals("00:00:00") && show_type.equals("3D")){
           return 4;
       }
       else{
           return 5;
       }
 }
}
